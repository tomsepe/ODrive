#include <math.h>
#include <algorithm>
#include <trapTraj.hpp>

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

// #include "odrive_main.h"
// #include "utils.h"

#define SQ(x) ((x) * (x))

// A sign function where input 0 has positive sign (not 0)
float sign_hard(float val)
{
    return (std::signbit(val)) ? -1.0f : 1.0f;
}

// Symbol                     Description
// Ta, Tv and Td              Duration of the stages of the AL profile
// Xi and Vi                  Adapted initial conditions for the AL profile
// Xf                         Position set-point
// s                          Direction (sign) of the trajectory
// Vmax, Amax, Dmax and jmax  Kinematic bounds
// Ar, Dr and Vr              Reached values of acceleration and velocity

TrapezoidalTrajectory::TrapezoidalTrajectory(Config_t &config) : config_(config) {}

bool TrapezoidalTrajectory::planTrapezoidal(float Xf, float Xi, float Vi,
                                            float Vmax, float Amax, float Dmax)
{
    float dX = Xf - Xi;                          // Distance to travel
    float stop_dist = (Vi * Vi) / (2.0f * Dmax); // Minimum stopping distance
    float dXstop = std::copysign(stop_dist, Vi); // Minimum stopping displacement
    float s = sign_hard(dX - dXstop);            // Sign of coast velocity (if any)
    Ar_ = s * Amax;                              // Maximum Acceleration (signed)
    Dr_ = -s * Dmax;                             // Maximum Deceleration (signed)
    Vr_ = s * Vmax;                              // Maximum Velocity (signed)

    // If we start with a speed faster than cruising, then we need to decel instead of accel
    // aka "double deceleration move" in the paper
    if ((s * Vi) > (s * Vr_))
    {
        // cout << "Handbrake!" << endl;
        Ar_ = -s * Amax;
    }

    // Time to accel/decel to/from Vr (cruise speed)
    Ta_ = (Vr_ - Vi) / Ar_;
    Td_ = -Vr_ / Dr_;

    // Integral of velocity ramps over the full accel and decel times to get
    // minimum displacement required to reach cuising speed
    float dXmin = 0.5f * Ta_ * (Vr_ + Vi) + 0.5f * Td_ * Vr_;

    // Are we displacing enough to reach cruising speed?
    if (s * dX < s * dXmin)
    {
        // cout << "Short Move: " << endl;
        // Short move (triangle profile)
        Vr_ = s * sqrtf((Dr_ * SQ(Vi) + 2 * Ar_ * Dr_ * dX) / (Dr_ - Ar_));
        Ta_ = std::max(0.0f, (Vr_ - Vi) / Ar_);
        Td_ = std::max(0.0f, -Vr_ / Dr_);
        Tv_ = 0.0f;
    }
    else
    {
        // cout << "Long Move: " << endl;
        // Long move (trapezoidal profile)
        Tv_ = (dX - dXmin) / Vr_;
    }

    // Fill in the rest of the values used at evaluation-time
    Tf_ = Ta_ + Tv_ + Td_;
    Xi_ = Xi;
    Xf_ = Xf;
    Vi_ = Vi;
    yAccel_ = Xi + Vi * Ta_ + 0.5f * Ar_ * SQ(Ta_); // pos at end of accel phase

    // cout << "Xi: " << Xi << "\tXf: " << Xf << "\tVi: " << Vi << endl;
    // cout << "Amax: " << Amax << "\tVmax: " << Vmax << "\tDmax: " << Dmax << endl;
    // cout << "dX: " << dX << "\tdXst: " << dXstop << "\tdXmin: " << dXmin << endl;
    // cout << "Ar: " << Ar_ << "\tVr: " << Vr_ << "\tDr: " << Dr_ << endl;
    // cout << "Ta: " << Ta_ << "\tTv: " << Tv_ << "\tTd: " << Td_ << endl;

    return true;
}

TrapezoidalTrajectory::Step_t TrapezoidalTrajectory::eval(float t)
{
    Step_t trajStep;
    if (t < 0.0f)
    { // Initial Condition
        trajStep.Y = Xi_;
        trajStep.Yd = Vi_;
        trajStep.Ydd = 0.0f;
    }
    else if (t < Ta_)
    { // Accelerating
        trajStep.Y = Xi_ + Vi_ * t + 0.5f * Ar_ * SQ(t);
        trajStep.Yd = Vi_ + Ar_ * t;
        trajStep.Ydd = Ar_;
    }
    else if (t < Ta_ + Tv_)
    { // Coasting
        trajStep.Y = yAccel_ + Vr_ * (t - Ta_);
        trajStep.Yd = Vr_;
        trajStep.Ydd = 0.0f;
    }
    else if (t < Tf_)
    { // Deceleration
        float td = t - Tf_;
        trajStep.Y = Xf_ + 0.5f * Dr_ * SQ(td);
        trajStep.Yd = Dr_ * td;
        trajStep.Ydd = Dr_;
    }
    else if (t >= Tf_)
    { // Final Condition
        trajStep.Y = Xf_;
        trajStep.Yd = 0.0f;
        trajStep.Ydd = 0.0f;
    }
    else
    {
        // TODO: report error here
    }

    return trajStep;
}

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <random>
#include <doctest.h>

using doctest::Approx;

// TEST_CASE("[planTrap]")
// {
//     TrapezoidalTrajectory::Config_t config;
//     TrapezoidalTrajectory trapTraj(config);

//     float pos_range = 10000.0;
//     float Vmax_range = 8000.0;
//     float Amax_range = 10000.0;

//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_real_distribution<float> posRandom(-pos_range, pos_range);
//     std::uniform_real_distribution<float> velRandom(0.1 * Vmax_range, Vmax_range);
//     std::uniform_real_distribution<float> accelRandom(0.1 * Amax_range, Amax_range);

//     std::uniform_real_distribution<float> viChance(0, 1);

//     cout << std::fixed;
//     cout << std::setprecision(2);
    
//     for (int i = 0; i < 10000; i++)
//     {
//         auto Vmax = velRandom(gen);
//         auto Amax = accelRandom(gen);
//         auto Dmax = Amax;

//         auto Xf = posRandom(gen);
//         auto Xi = posRandom(gen);

//         float Vi = 0.0f;
//         std::uniform_real_distribution<float> viRandom(-Vmax * 1.5, Vmax * 1.5);
//         if(viChance(gen) <= 0.5)
//             Vi = viRandom(gen);
//         else
//             Vi = 0.0f;

//         cout << endl;
//         cout << "Test " << i+1 << endl;
//         trapTraj.planTrapezoidal(Xf, Xi, Vi, Vmax, Amax, Dmax);
        
//         TrapezoidalTrajectory::Step_t step;

//         step = trapTraj.eval(0.0f);
//         CHECK(step.Y == Approx(Xi).epsilon(0.0001));
//         CHECK(step.Yd == Approx(Vi).epsilon(0.0001));

//         step = trapTraj.eval(trapTraj.Tf_ - 0.001);
//         CHECK(step.Y == Approx(Xf-0.001*step.Yd).epsilon(0.01));
//         CHECK(step.Yd == Approx(0.0f-Dmax*0.001).epsilon(0.01));

//         step = trapTraj.eval(trapTraj.Tf_);
//         CHECK(step.Y == Approx(Xf).epsilon(0.0001));
//         CHECK(step.Yd == Approx(0.0f).epsilon(0.0001));
//     }



//     cout << "\n\n\n" << endl;
// }

TEST_CASE("handbrake"){
    TrapezoidalTrajectory::Config_t config;
    TrapezoidalTrajectory trapTraj(config);



    float Xi = 7310.54;
    float Xf = -9379.76;
    float Vi = 3319.07;

    config.accel_limit = 3711.41;
    config.decel_limit = config.accel_limit;
    config.vel_limit = 5055.92;

    trapTraj.planTrapezoidal(Xf, Xi, Vi, config.vel_limit, config.accel_limit, config.decel_limit);
    float dT = trapTraj.Tf_ / 1000.0f;
    float time = 0;

    float lastTime = 0.0f;
    for(int i = 0; i < 1000; ++i){
        auto step = trapTraj.eval(time);
        cout << time+lastTime << "," << step.Y << "," << step.Yd << "," << step.Ydd << "\n";
        if(i == 500){
            lastTime = time;
            Xi = step.Y;
            Xf = 15000.0;
            Vi = step.Yd;
            trapTraj.planTrapezoidal(Xf, Xi, Vi, config.vel_limit, config.accel_limit, config.decel_limit);
            dT = trapTraj.Tf_ / 500.0f;
            time = 0.0f;
        }
        time += dT;
    }
}