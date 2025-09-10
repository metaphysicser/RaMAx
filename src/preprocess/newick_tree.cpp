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

// ⚠️ 工具：返回两节点距离（用父指针一路向上，适合叶子数不大）
double NewickParser::distanceBetween(int u, int v) const {
    // 记录到根的路径及累积距离
    std::unordered_map<int, double> distUp;
    double acc = 0.0;
    int x = u;
    while (x != -1) {
        distUp[x] = acc;
        int p = nodes_[x].father;
        if (p == -1) break;
        acc += nodes_[x].branchLength;
        x = p;
    }
    // 再从 v 向上找最近公共祖先
    double accV = 0.0;
    int y = v;
    while (y != -1) {
        if (distUp.count(y)) {
            return accV + distUp[y];
        }
        int p = nodes_[y].father;
        if (p == -1) break;
        accV += nodes_[y].branchLength;
        y = p;
    }
    return acc + accV; // 保险：不应走到这里
}

// --------------------------------------------------
// 迭代贪心中心排序
// --------------------------------------------------
std::vector<int> NewickParser::orderLeavesGreedyMinSum(int rootLeaf) {
    // -------- 1. 收集所有叶节点 --------
    std::vector<int> leaves;
    leaves.reserve(nodes_.size());
    for (const auto& n : nodes_) if (n.isLeaf) leaves.push_back(n.id);
    const int L = static_cast<int>(leaves.size());
    if (L == 0) return {};

    // -------- 2. 预计算 L×L 距离矩阵 --------
    std::vector<std::vector<double>> D(L, std::vector<double>(L, 0.0));
    for (int i = 0; i < L; ++i) {
        for (int j = i + 1; j < L; ++j) {
            double d = distanceBetween(leaves[i], leaves[j]);
            D[i][j] = D[j][i] = d;
        }
    }

    // 映射 leafID -> 索引
    std::unordered_map<int, int> leaf2idx;
    for (int i = 0;i < L;++i) leaf2idx[leaves[i]] = i;

    // -------- 3. 初始化 sumDist --------
    std::vector<double> sumDist(L, 0.0);
    for (int i = 0;i < L;++i)
        for (int j = 0;j < L;++j) if (i != j) sumDist[i] += D[i][j];

    // -------- 4. 迭代贪心选择 --------
    std::vector<int> order; order.reserve(L);
    std::vector<char> removed(L, 0);   // 0=在集合，1=已移除

    for (int step = 0; step < L; ++step) {
        // 找剩余里 sumDist 最小者
        int bestIdx = -1;
        double bestVal = std::numeric_limits<double>::max();
        for (int i = 0;i < L;++i) {
            if (removed[i]) continue;
            if (sumDist[i] < bestVal) {
                bestVal = sumDist[i];
                bestIdx = i;
            }
        }
        // 记录并移除
        order.push_back(leaves[bestIdx]);
        removed[bestIdx] = 1;

        // 更新其余节点的 sumDist
        for (int j = 0;j < L;++j) {
            if (removed[j]) continue;
            sumDist[j] -= D[j][bestIdx];
        }
    }
    return order;  // 下标 0 是全局距离最小中心，后面依次次优
}

std::vector<std::string> NewickParser::getLeafNames() const {
    std::vector<std::string> leaves;
    leaves.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        if (node.isLeaf) {
            // 只有名称非空才加入
            if (!node.name.empty()) {
                leaves.push_back(node.name);
            }
        }
    }
    return leaves;
}

// 在 NewickParser 声明里加（可放 public）：
int NewickParser::findNodeIdByName(const std::string& name) const {
    for (const auto& n : nodes_) {
        if (n.name == name) return n.id;
    }
    return -1;
}

// 传入根 id，收集该子树所有 oldId，重建 nodes_ 并重新编号
void NewickParser::restrictToSubtreeByRootId(int rootId) {
    if (rootId < 0 || rootId >= (int)nodes_.size())
        throw std::runtime_error("restrictToSubtreeByRootId: invalid root id");

    // 收集子树所有 oldId
    std::vector<int> stack{ rootId };
    std::vector<int> kept; kept.reserve(nodes_.size());
    std::vector<char> mark(nodes_.size(), 0);

    while (!stack.empty()) {
        int u = stack.back(); stack.pop_back();
        if (mark[u]) continue;
        mark[u] = 1;
        kept.push_back(u);
        const auto& nd = nodes_[u];
        if (nd.leftChild != -1) stack.push_back(nd.leftChild);
        if (nd.rightChild != -1) stack.push_back(nd.rightChild);
    }

    // oldId -> newId 映射
    std::vector<int> old2new(nodes_.size(), -1);
    std::vector<NewickTreeNode> newNodes; newNodes.reserve(kept.size());

    // 先按 kept 的顺序建立新节点（重新编号为 0..k-1）
    for (size_t i = 0; i < kept.size(); ++i) {
        int oldId = kept[i];
        old2new[oldId] = static_cast<int>(i);
        NewickTreeNode nn = nodes_[oldId];
        nn.id = static_cast<int>(i);
        newNodes.push_back(std::move(nn));
    }

    // 修复 parent/children 指针到新编号；把 root 的 father 设为 -1
    for (auto& nn : newNodes) {
        // 映射 child
        nn.leftChild = (nn.leftChild != -1 && old2new[nn.leftChild] != -1) ? old2new[nn.leftChild] : -1;
        nn.rightChild = (nn.rightChild != -1 && old2new[nn.rightChild] != -1) ? old2new[nn.rightChild] : -1;
        // 映射 father
        nn.father = (nn.father != -1 && old2new[nn.father] != -1) ? old2new[nn.father] : -1;
    }

    // 确保新 root 的 father == -1
    int newRoot = old2new[rootId];
    if (newRoot >= 0 && newRoot < (int)newNodes.size()) {
        newNodes[newRoot].father = -1;
    }

    nodes_.swap(newNodes);
    currentIndex_ = static_cast<int>(nodes_.size());
}
