#include "data_process.h"

#include <vector>
#include <stdexcept>
#include <cctype>    // std::isspace
#include <cstdlib>   // std::strtod
#include <cstring>   // std::strlen
#include <sstream>   // (for future extensions if needed)

NewickParser::NewickParser(const std::string& newickStr) {
    parse(newickStr);
}

const std::vector<NewickTreeNode>& NewickParser::getNodes() const {
    return nodes_;
}

void NewickParser::parse(const std::string& newickStr) {
    // Copy the string into a modifiable character array for pointer-based parsing.
    std::vector<char> buffer(newickStr.begin(), newickStr.end());
    buffer.push_back('\0'); // Null-terminate.

    int index = 0;
    int length = static_cast<int>(buffer.size()) - 1; // Exclude null terminator.

    // Start parsing from the root with father = -1
    int rootId = parseSubtree(buffer.data(), index, length, -1);

    // After parsing the subtree, skip whitespace and expect a semicolon ';'
    skipWhitespace(buffer.data(), index, length);
    if (index < length && buffer[index] == ';') {
        index++;
        skipWhitespace(buffer.data(), index, length);
    }

    // If there are remaining non-whitespace characters, the format is invalid
    if (index < length) {
        throw std::runtime_error("Newick format error: Extra characters after parsing.");
    }

    // If the root node has no name, assign a default name "root"
    if (rootId >= 0 && nodes_[rootId].name.empty()) {
        nodes_[rootId].name = "root";
    }
}

int NewickParser::parseSubtree(char* str, int& index, int length, int father) {
    skipWhitespace(str, index, length);

    // Create a new node
    NewickTreeNode node;
    node.id = currentIndex_++;
    node.father = father;
    node.branchLength = 0.0;
    node.isLeaf = false;
    node.leftChild = -1;
    node.rightChild = -1;

    // Append a placeholder; we'll overwrite fields as we parse children/leaf.
    nodes_.push_back(node);

    skipWhitespace(str, index, length);

    // Check if this is an internal node (starts with '(')
    if (index < length && str[index] == '(') {
        // Skip '('
        index++;
        skipWhitespace(str, index, length);

        // Parse left child
        int leftId = parseSubtree(str, index, length, node.id);
        nodes_[node.id].leftChild = leftId;

        skipWhitespace(str, index, length);

        // Expect a comma
        if (index >= length || str[index] != ',') {
            throw std::runtime_error("Newick format error: Expected ',' after left subtree.");
        }
        index++;  // Skip comma
        skipWhitespace(str, index, length);

        // Parse right child
        int rightId = parseSubtree(str, index, length, node.id);
        nodes_[node.id].rightChild = rightId;

        skipWhitespace(str, index, length);

        // Expect closing ')'
        if (index >= length || str[index] != ')') {
            throw std::runtime_error("Newick format error: Expected ')' after right subtree.");
        }
        index++;  // Skip ')'
        skipWhitespace(str, index, length);

        // Parse optional internal node name
        parseNodeName(str, index, length, nodes_[node.id].name);

        // Parse optional branch length for this internal node
        if (index < length && str[index] == ':') {
            index++;
            skipWhitespace(str, index, length);
            nodes_[node.id].branchLength = parseBranchLength(str, index, length);
        }

        // It's confirmed as an internal node
        nodes_[node.id].isLeaf = false;
    }
    else {
        // It's a leaf node
        nodes_[node.id].isLeaf = true;

        // Parse leaf name
        parseNodeName(str, index, length, nodes_[node.id].name);

        // Parse optional branch length
        if (index < length && str[index] == ':') {
            index++;
            skipWhitespace(str, index, length);
            nodes_[node.id].branchLength = parseBranchLength(str, index, length);
        }

        nodes_[node.id].isLeaf = true;
    }

    return node.id;
}

void NewickParser::parseNodeName(char* str, int& index, int length, std::string& outName) {
    skipWhitespace(str, index, length);

    // Node names can include alphanumeric characters and underscores,
    // but stop at ':', ',', ')', '(', ';', or whitespace.
    while (index < length) {
        char c = str[index];
        if (c == ':' || c == ',' || c == ')' || c == '(' || c == ';' || std::isspace(static_cast<unsigned char>(c))) {
            break;
        }
        outName.push_back(c);
        index++;
    }

    // Trim any leading/trailing whitespace just in case
    trimString(outName);
}

double NewickParser::parseBranchLength(char* str, int& index, int length) {
    skipWhitespace(str, index, length);
    char* startPtr = str + index;
    char* endPtr = nullptr;

    // Convert substring to double
    double val = std::strtod(startPtr, &endPtr);
    if (startPtr == endPtr) {
        throw std::runtime_error("Newick format error: Invalid branch length.");
    }

    // Advance index by number of characters consumed
    int consumed = static_cast<int>(endPtr - startPtr);
    index += consumed;

    skipWhitespace(str, index, length);
    return val;
}

void NewickParser::skipWhitespace(char* str, int& index, int length) {
    while (index < length && std::isspace(static_cast<unsigned char>(str[index]))) {
        index++;
    }
}

void NewickParser::trimString(std::string& s) {
    if (s.empty()) return;

    // Trim leading whitespace
    std::size_t startPos = 0;
    while (startPos < s.size() && std::isspace(static_cast<unsigned char>(s[startPos]))) {
        startPos++;
    }

    // Trim trailing whitespace
    std::size_t endPos = s.size();
    while (endPos > startPos && std::isspace(static_cast<unsigned char>(s[endPos - 1]))) {
        endPos--;
    }

    s = s.substr(startPos, endPos - startPos);
}
