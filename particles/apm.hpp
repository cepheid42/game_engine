#ifndef GAME_ENGINE_APM_HPP
#define GAME_ENGINE_APM_HPP

#include "particles.hpp"

#include <cassert>
#include <map>
#include <span>
#include <vector>

namespace tf::apm {

template<int DIM=3>
class KdTree {
   struct Node {
      std::array<double, DIM> point;
      Node* left{nullptr};
      Node* right{nullptr};
   };

   Node* root{nullptr};

   auto insertRecursive(Node* node, const auto& point, const int depth) -> Node* {
      if (node == nullptr) { return new Node{point}; }

      const auto cd = depth % DIM;

      if (point[cd] < node->point[cd]) {
         node->left = insert(node->left, point, depth + 1);
      } else {
         node->right = insert(node->right, point, depth + 1);
      }
      return node;
   }

   auto searchRecursive(Node* node, const auto& point, const int depth) const -> bool {
      if (node == nullptr) { return false; }
      if (node->point == point) { return true; } // std::arrays don't have comparisons do they?

      const auto cd = depth % DIM;

      if (point[cd] < node->point[cd]) {
         return search(node->left, point, depth + 1);
      } else {
         return search(node->right, point, depth + 1);
      }
   }

   void freeTree(Node* node) {
      if (!node) { return; }
      freeTree(node->left);
      freeTree(node->right);
      delete node;
   }

   void printRecursive(Node* node, const int depth) const {
      if (node == nullptr) { return; }
      // Print current node with indentation based on depth
      for (auto i = 0; i < depth; i++) std::cout << "  ";
      std::cout << "(";
      for (size_t i = 0; i < DIM; i++) {
         std::cout << node->point[i];
         if (i < DIM - 1) std::cout << ", ";
      }
      std::cout << ")" << std::endl;

      // Recursively print left and right children
      printRecursive(node->left, depth + 1);
      printRecursive(node->right, depth + 1);
   }

   void nearestNeighbor(Node* node, const int depth, auto best) {
      if (node == nullptr) { return; }
      // if (())


   }

public:
   KdTree() = default;

   ~KdTree() {
      freeTree(root);
   }

   void insert(const auto& point) { root = insertRecursive(root, point, 0); }
   auto search(const auto& point) -> bool { return searchRecursive(root, point, 0); }
   void print() const { printRecursive(root, 0); }
};


} // end namespace tf::apm
#endif //GAME_ENGINE_APM_HPP
