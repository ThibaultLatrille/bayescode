#include "Chronogram.hpp"
#include <cassert>
#include <fstream>
#include <sstream>

NodeAges::NodeAges(const Tree &intree, const string &fossilsfile)
    : SimpleNodeArray<double>(intree, 0.0) {
    if (fossilsfile != "Null") {
        // for line in file
        std::ifstream input_stream{fossilsfile};
        if (!input_stream) {
            std::cerr << "Fossils file " << fossilsfile << " doesn't exist" << std::endl;
        }

        std::string line;
        // skip the header of the file
        getline(input_stream, line);
        char sep{'\t'};
        {
            std::istringstream line_stream(line);
            std::string word;
            // Skip the first column (taxon name)
            getline(line_stream, word, sep);
            assert(word == "NodeName");
            getline(line_stream, word, sep);
            assert(word == "Age");
            getline(line_stream, word, sep);
            assert(word == "LowerBound");
            getline(line_stream, word, sep);
            assert(word == "UpperBound");
        }
        while (getline(input_stream, line)) {
            std::istringstream line_stream(line);
            std::string word;
            // Skip the first column (taxon name)
            Tree::NodeIndex node_clamped{-1};
            getline(line_stream, word, sep);
            if (word == "Root") {
                node_clamped = GetTree().root();
            } else {
                for (auto const &node : GetTree().root_to_leaves_iter()) {
                    if (GetTree().node_name(node) == word) { node_clamped = node; }
                }
            }
            if (node_clamped == -1) {
                std::cerr << word << " node is unknown." << std::endl;
                continue;
            }
            node_clamped_set.insert(node_clamped);
            getline(line_stream, word, sep);
            clamped_ages[node_clamped] = std::stod(word);
            (*this)[node_clamped] = clamped_ages[node_clamped];
            getline(line_stream, word, sep);
            clamped_lower_bound[node_clamped] = std::stod(word);
            assert(clamped_lower_bound[node_clamped] <= clamped_ages[node_clamped]);
            getline(line_stream, word, sep);
            clamped_upper_bound[node_clamped] = std::stod(word);
            assert(clamped_upper_bound[node_clamped] >= clamped_ages[node_clamped]);
        }
        if (node_clamped_set.count(GetTree().root()) == 0) {
            double max_root_age = 0.0;
            double root_ecc = EccentricityRecursive(GetTree().root());
            for (Tree::NodeIndex const &node_clamped : node_clamped_set) {
                double d =
                    this->GetVal(node_clamped) * root_ecc / EccentricityRecursive(node_clamped);
                if (d > max_root_age) { max_root_age = d; }
            }
            (*this)[GetTree().root()] = max_root_age;
        }

    } else {
        (*this)[GetTree().root()] = 1.0;
        node_clamped_set.insert(GetTree().root());
        clamped_ages[GetTree().root()] = 1.0;
        clamped_lower_bound[GetTree().root()] = 1.0;
        clamped_upper_bound[GetTree().root()] = 1.0;
    }
    std::set<Tree::NodeIndex> node_init_set = node_clamped_set;
    if (node_clamped_set.count(GetTree().root()) == 0) { node_init_set.insert(GetTree().root()); }
    for (Tree::NodeIndex const &node_init : node_init_set) {
        std::set<Tree::NodeIndex> internal_node_set{};
        std::set<Tree::NodeIndex> frontier_node_set{};
        std::unordered_map<Tree::NodeIndex, double> distance{{node_init, 0.0}};

        std::queue<Tree::NodeIndex> bfs_queue{};
        bfs_queue.push(node_init);
        while (!bfs_queue.empty()) {
            Tree::NodeIndex node_bfs = bfs_queue.front();
            bfs_queue.pop();
            if ((node_clamped_set.count(node_bfs) == 1 and node_bfs != node_init) or
                GetTree().is_leaf(node_bfs)) {
                frontier_node_set.insert(node_bfs);
            } else {
                if (node_bfs != node_init) { internal_node_set.insert(node_bfs); }
                for (Tree::NodeIndex child : GetTree().children(node_bfs)) {
                    distance[child] = distance[node_bfs] + 1;
                    bfs_queue.push(child);
                }
            }
        }

        // Update the partial subtree (until the clamped frontier)
        double min_delta_t = std::numeric_limits<double>::infinity();
        for (Tree::NodeIndex const &frontier_node : frontier_node_set) {
            double frontier_age =
                node_clamped_set.count(frontier_node) == 1 ? clamped_ages[frontier_node] : 0.0;
            double delta_t = (clamped_ages[node_init] - frontier_age) / distance[frontier_node];
            if (delta_t < min_delta_t) { min_delta_t = delta_t; }
        }

        for (Tree::NodeIndex const &internal_node : internal_node_set) {
            assert((*this)[GetTree().parent(internal_node)] != 0.0);
            (*this)[internal_node] = (*this)[GetTree().parent(internal_node)] - min_delta_t;
            assert((*this)[internal_node] > 0.0);
        }
    }
    assert(Check());
}

double NodeAges::EccentricityRecursive(Tree::NodeIndex node) {
    if (GetTree().is_leaf(node)) {
        // set leaf nodes at age 0
        (*this)[node] = 0.0;
    } else {
        // proceed from leaves to root using recursive algorithm
        double max_eccent = 0.0;
        for (auto const &child : GetTree().children(node)) {
            double eccent = EccentricityRecursive(child) + 1.0;
            if (eccent > max_eccent) { max_eccent = eccent; }
        }

        (*this)[node] = max_eccent;
    }
    return GetVal(node);
}

bool NodeAges::Check() const {
    for (auto const &node : GetTree().root_to_leaves_iter()) {
        if (!GetTree().is_root(node)) {
            if (GetTree().is_leaf(node)) {
                assert(GetVal(node) == 0.0);
            } else {
                assert(GetVal(node) != 0.0);
            }
            assert(GetVal(GetTree().parent(node)) >= GetVal(node));
        }
    }
    return true;
}

void NodeAges::SlidingMove(Tree::NodeIndex node, double scale) {
    assert(!GetTree().is_leaf(node));

    double lower_bound = node_clamped_set.count(node) == 1 ? clamped_lower_bound[node] : 0.0;
    for (auto const &child : GetTree().children(node)) {
        if (GetVal(child) > lower_bound) { lower_bound = GetVal(child); }
    }
    double upper_bound = node_clamped_set.count(node) == 1 ? clamped_upper_bound[node] : 1.0;
    if (!GetTree().is_root(node)) {
        double p = GetVal(GetTree().parent(node));
        if (p < upper_bound) { upper_bound = p; }
    }

    assert(upper_bound >= GetVal(node));
    assert(GetVal(node) >= lower_bound);
    if (upper_bound == lower_bound) { return; }

    double x = GetVal(node) + scale * (upper_bound - lower_bound);

    while ((x < lower_bound) || (x > upper_bound)) {
        if (x < lower_bound) { x = 2 * lower_bound - x; }
        if (x > upper_bound) { x = 2 * upper_bound - x; }
    }

    (*this)[node] = x;
}

Chronogram::Chronogram(const NodeAges &innodeages)
    : SimpleBranchArray<double>(innodeages.GetTree()), nodeages(innodeages) {
    Update();
}

void Chronogram::Update() {
    for (auto const &node : GetTree().root_to_leaves_iter()) {
        if (!GetTree().is_root(node)) { UpdateBranch(GetTree().parent(node), node); }
    }
}

void Chronogram::UpdateLocal(Tree::NodeIndex node) {
    // update all branch lengths around this node

    // for the branch attached to the node
    if (GetTree().is_root(node)) {
        Update();
    } else {
        UpdateBranch(GetTree().parent(node), node);
        // for all children
        for (auto const &child : GetTree().children(node)) { UpdateBranch(node, child); }
    }
}

void Chronogram::UpdateBranch(Tree::NodeIndex parent, Tree::NodeIndex node) {
    (*this)[GetTree().branch_index(node)] =
        (nodeages.GetVal(parent) - nodeages.GetVal(node)) / nodeages.GetVal(GetTree().root());
}
