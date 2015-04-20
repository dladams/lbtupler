// gettrees.h

#include <string>

class TTree;

// Following functions return pointers to trees.
// If provided, the name is used to locate the tree.
// Use tname = "null" to zero the pointer.
TTree* simtree(std::string tname ="");
