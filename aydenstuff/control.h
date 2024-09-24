//
// Created by akis on 8/26/24.
//

#ifndef TRIFORCE_CONTROL_H
#define TRIFORCE_CONTROL_H

#include <iostream>
#include <variant>
#include <concepts>
#include <source_location>

#include "error_functions.h"

namespace control {
  enum class IntervalType : unsigned { WallTime=0, SimTime=1, StepCount=2 };
  
  inline std::ostream& operator<<(std::ostream& os, const IntervalType& type)
  {
    auto output_string = std::string();
    
    switch (type) {
      case IntervalType::WallTime:
        output_string = "WallTime";
        break;
      case IntervalType::SimTime:
        output_string = "SimTime";
        break;
      case IntervalType::StepCount:
        output_string = "StepCount";
        break;
    }
    
    os << output_string;
    return os;
    //
  } // end operator<<(std::ostream&, const IntervalType&)
  
  template <std::floating_point fp>
  struct Interval {
    IntervalType type;
    std::variant<size_t, fp> value;
    std::variant<size_t, fp> counter;
  };
  
  template <std::floating_point fp>
  inline std::ostream& operator<<(std::ostream& os, const Interval<fp>& interval)
  {
    os << "    type: " << interval.type << "\n";
    std::visit([&](auto&& arg) {os << "    value: " << arg;}, interval.value);
    
    return os;
  }
  
  template<std::floating_point fp, std::integral I>
  bool updateThisStep(Interval<fp> &interval, [[maybe_unused]] I step) {
    if (step == 0) { return true; }
    
    // Currently, if T is size_t, then interval.type == StepCount
    const size_t value = std::get<I>(interval.value);
    auto& counter = std::get<I>(interval.counter);
    
    counter++; // <AJK> A little spooky, but we'll give it a shot
    if (counter >= value) {
      counter = 0u;
      return true;
    }
    
    return false;
    //
  } // end updateThisStep<size_t>
  
  template<std::floating_point fp>
  bool updateThisStep(Interval<fp> &interval, fp step) {
    if (step == 0) { return true; }
    
    const auto& value = std::get<fp>(interval.value);
    auto& counter = std::get<fp>(interval.counter);
    
    switch (interval.type) {
      case IntervalType::StepCount:
        ERROR("Mismatched types in interval query.", std::source_location::current());
        break;
      case IntervalType::WallTime:
        ERROR("Wall time interval not yet implemented.", std::source_location::current());
        break;
      case IntervalType::SimTime:
        counter += step;
        if (counter >= value) {
          counter = static_cast<fp>(0.0);
          return true;
        } else {
          return false;
        }
        break;
    }
    
    return false;
    //
  } // end updateThisStep<fptype>
  //
} // end namespace control
#endif //TRIFORCE_CONTROL_H
