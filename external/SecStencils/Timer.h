#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <vector>
#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "color.h"

class Timer
{
public:
    static
    void Start(
        std::vector<double>& timer,
        const int color,
        const std::string msg)
    {
        SetColor(color);
        std::cout << msg << white << " ... " << "\n";
        #ifdef _OPENMP
            timer.push_back(omp_get_wtime());
        #else
            timer.push_back(clock());
        #endif
    }

    static
    void Start(
        std::vector<double>& timer,
        const int color,
        const std::string msg,
        const double val)
    {
        SetColor(color);
        std::cout << msg << white;
        std::cout << " (" << val << ") ... " << "\n";
        #ifdef _OPENMP
            timer.push_back(omp_get_wtime());
        #else
            timer.push_back(clock());
        #endif
    }

    static
    void Stop(
        std::vector<double>& timer,
        const int color,
        const std::string msg = std::string())
    {
        double duration = GetDuration(timer.back());
        timer.pop_back();
        SetColor(color);
        std::cout 
        << "done" << white
        << " (" << duration << " s) " 
        << msg << "\n";
    }

    static
    void SetColor(const int color)
    {
        switch (color)
        {
            case COLOR_WHITE:
                std::cout << white;
                break;
            case COLOR_RED:
                std::cout << red;
                break;
            case COLOR_BLUE:
                std::cout << blue;
                break;
            case COLOR_GREEN:
                std::cout << green;
                break;
            case COLOR_YELLOW:
                std::cout << yellow;
                break;
        }
    }

    static
    double 
    GetDuration(const double init)
    {
        #ifdef _OPENMP
            return (omp_get_wtime() - init);
        #else
            return (clock() - init) / CLOCKS_PER_SEC;
        #endif
    }
};

#endif // TIMER_H
