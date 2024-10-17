//
// Created by cepheid on 9/27/24.
//

#ifndef BREAD_H
#define BREAD_H

struct DefaultPolicy1 {};
struct DefaultPolicy2 {};
struct DefaultPolicy3 {};
struct DefaultPolicy4 {};

struct CustomPolicy {
  static void doPrint() { std::cout << "CustomPolicy" << std::endl; }
};

struct DefaultPolicies {
  using P1 = DefaultPolicy1;
  using P2 = DefaultPolicy1;
  using P3 = DefaultPolicy1;
  using P4 = DefaultPolicy1;
};

template<typename Base, int D>
struct Discriminator : Base {};

template<typename Setter1, typename Setter2, typename Setter3, typename Setter4>
struct PolicySelector : Discriminator<Setter1, 1>,
                        Discriminator<Setter2, 2>,
                        Discriminator<Setter3, 3>,
                        Discriminator<Setter4, 4>
{};

struct DefaultPolicyArgs : virtual DefaultPolicies {};

template<typename Policy>
struct Policy1_is : virtual DefaultPolicies {
  using P1 = Policy;
};

template<typename Policy>
struct Policy2_is : virtual DefaultPolicies {
  using P2 = Policy;
};

template<typename Policy>
struct Policy3_is : virtual DefaultPolicies {
  using P3 = Policy;
};

template<typename Policy>
struct Policy4_is : virtual DefaultPolicies {
  using P4 = Policy;
};

template<typename PolicySetter1 = DefaultPolicyArgs,
         typename PolicySetter2 = DefaultPolicyArgs,
         typename PolicySetter3 = DefaultPolicyArgs,
         typename PolicySetter4 = DefaultPolicyArgs>
struct BreadSlicer {
  using Policies = PolicySelector<PolicySetter1, PolicySetter2, PolicySetter3, PolicySetter4>;
  using P3 = typename Policies::P3;

  void print() {
    P3::doPrint();
  }
};



#endif //BREAD_H
