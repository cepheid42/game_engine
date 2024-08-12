//
// Created by cepheid on 7/31/24.
//

#ifndef PROFILING_H
#define PROFILING_H


// struct AutoProfile
// {
//   AutoProfile(const char* name)
//   {
//     m_name = name;
//     m_startTime = QueryPerformanceCounter();
//   }
//
//   ~AutoProfile()
//   {
//     std::int64_t endTime = QueryPerformanceCounter();
//     std::int64_t elapsedTime = endTime - m_startTime;
//     g_profileManager.storeSample(m_name, elapsedTime);
//   }
//   const char*
//   std::int64_t
//   m_name;
//   m_startTime;
// };
// #define PROFILE(name) AutoProfile p(name)


#endif //PROFILING_H
