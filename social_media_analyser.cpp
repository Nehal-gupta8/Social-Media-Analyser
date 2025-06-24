#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

class SocialGraph {
private:
    unordered_map<string, vector<string>> adjList;       // User -> Following
    unordered_map<string, vector<string>> incomingList; // User -> Followers
    unordered_set<string> users;

public:
    void addUser(const string& user) {
        if (users.find(user) == users.end()) {
            users.insert(user);
            adjList[user] = vector<string>();
            incomingList[user] = vector<string>();
        }
    }

    void addConnection(const string& from, const string& to) {
        addUser(from);
        addUser(to);
        
        // Add to following list
        if (find(adjList[from].begin(), adjList[from].end(), to) == adjList[from].end()) {
            adjList[from].push_back(to);
        }
        
        // Add to followers list
        if (find(incomingList[to].begin(), incomingList[to].end(), from) == incomingList[to].end()) {
            incomingList[to].push_back(from);
        }
    }

    const vector<string>& getFollowing(const string& user) const {
        static const vector<string> empty;
        auto it = adjList.find(user);
        return it != adjList.end() ? it->second : empty;
    }

    const vector<string>& getFollowers(const string& user) const {
        static const vector<string> empty;
        auto it = incomingList.find(user);
        return it != incomingList.end() ? it->second : empty;
    }

    const unordered_set<string>& getAllUsers() const {
        return users;
    }

    size_t getUserCount() const {
        return users.size();
    }

    const unordered_map<string, vector<string>>& getAdjacencyList() const {
        return adjList;
    }
};

class CommunityDetector {
private:
    const SocialGraph& graph;
    unordered_map<string, int> communities;
    double totalWeight;

public:
    CommunityDetector(const SocialGraph& g) : graph(g) {
        totalWeight = 0.0;
        // Initialize each node in its own community
        int communityId = 0;
        for (const auto& user_entry : graph.getAllUsers()) {
            communities[user_entry] = communityId++;
        }
    }

    // Calculate modularity gain if node is moved to a community
    double calculateModularityGain(const string& node, int newCommunity) {
        double gain = 0.0;
        int oldCommunity = communities[node];
        
        // Get all neighbors
        const auto& following = graph.getFollowing(node);
        const auto& followers = graph.getFollowers(node);
        
        unordered_set<string> neighbors;
        for (const auto& n : following) neighbors.insert(n);
        for (const auto& n : followers) neighbors.insert(n);
        
        // Calculate sum of weights from node to new community and old community
        double sumInNew = 0.0, sumInOld = 0.0;
        double kin = following.size() + followers.size(); // total degree of node
        
        for (const auto& neighbor : neighbors) {
            if (communities[neighbor] == newCommunity) {
                sumInNew += 1.0; // unweighted graph
            }
            if (communities[neighbor] == oldCommunity) {
                sumInOld += 1.0;
            }
        }
        
        // Calculate sum of weights in new and old communities
        double sumTotNew = 0.0, sumTotOld = 0.0;
        for (const auto& comm_entry : communities) {
            const string& user = comm_entry.first;
            int comm = comm_entry.second;
            if (comm == newCommunity) {
                sumTotNew += graph.getFollowing(user).size() + graph.getFollowers(user).size();
            }
            if (comm == oldCommunity) {
                sumTotOld += graph.getFollowing(user).size() + graph.getFollowers(user).size();
            }
        }
        
        // Modularity gain formula (simplified for directed graphs)
        gain = (sumInNew - sumInOld) / totalWeight - 
               kin * (sumTotNew - sumTotOld) / (2 * totalWeight * totalWeight);
        
        return gain;
    }

    unordered_map<int, vector<string>> detectCommunities(int maxIterations = 100) {
    bool improved = true;
    totalWeight = calculateTotalWeight();
    int iteration = 0;
    
    if (totalWeight == 0) {
        unordered_map<int, vector<string>> result;
        int i = 0;
        for (const auto& user : graph.getAllUsers()) {
            result[i++].push_back(user);
        }
        return result;
    }
    
    while (improved && iteration < maxIterations) {
        improved = false;
        iteration++;
        
        cout << "Community detection iteration " << iteration << endl; // Debug output
        
        // Phase 1: Modularity optimization
        for (const auto& node_entry : graph.getAllUsers()) {
            const string& node = node_entry;
            int bestCommunity = communities[node];
            double maxGain = 0.0;
            
            unordered_set<int> neighborCommunities;
            const auto& following = graph.getFollowing(node);
            const auto& followers = graph.getFollowers(node);
            
            for (const auto& n : following) neighborCommunities.insert(communities[n]);
            for (const auto& n : followers) neighborCommunities.insert(communities[n]);
            
            for (int comm : neighborCommunities) {
                if (comm != communities[node]) {
                    double gain = calculateModularityGain(node, comm);
                    if (gain > maxGain) {
                        maxGain = gain;
                        bestCommunity = comm;
                    }
                }
            }
            
            if (maxGain > 0.0 && bestCommunity != communities[node]) {
                communities[node] = bestCommunity;
                improved = true;
            }
        }
    }
    
    if (iteration >= maxIterations) {
        cout << "Warning: Reached maximum iterations (" << maxIterations << ") in community detection" << endl;
    }
    
    unordered_map<int, vector<string>> result;
    for (const auto& comm_entry : communities) {
        result[comm_entry.second].push_back(comm_entry.first);
    }
    
    return result;
}

private:
    double calculateTotalWeight() {
        double weight = 0.0;
        for (const auto& adj_entry : graph.getAdjacencyList()) {
            weight += adj_entry.second.size();
        }
        return weight;
    }
};

class SocialAnalyzer {
private:
    const SocialGraph& graph;

public:
    SocialAnalyzer(const SocialGraph& g) : graph(g) {}

    // PageRank for influencer detection
    unordered_map<string, double> calculatePageRank(double damping = 0.85, int iterations = 100, double tolerance = 1e-6) {
        unordered_map<string, double> pageRank;
        const size_t userCount = graph.getUserCount();
        if (userCount == 0) return pageRank;

        // Initialize ranks
        const double initialRank = 1.0 / userCount;
        for (const auto& user : graph.getAllUsers()) {
            pageRank[user] = initialRank;
        }

        for (int i = 0; i < iterations; ++i) {
            unordered_map<string, double> newRank;
            double totalDiff = 0.0;

            for (const auto& user_entry : graph.getAllUsers()) {
                const string& user = user_entry;
                double sum = 0.0;

                // Sum contributions from followers
                for (const auto& follower : graph.getFollowers(user)) {
                    const auto& following = graph.getFollowing(follower);
                    if (!following.empty()) {
                        sum += pageRank[follower] / following.size();
                    }
                }

                // Calculate new rank
                newRank[user] = (1.0 - damping) / userCount + damping * sum;
                totalDiff += abs(newRank[user] - pageRank[user]);
            }

            pageRank = newRank;

            // Check for convergence
            if (totalDiff < tolerance) {
                break;
            }
        }

        return pageRank;
    }

    // Jaccard similarity for friend recommendation
    double jaccardSimilarity(const string& user1, const string& user2) {
        const auto& following1 = graph.getFollowing(user1);
        const auto& following2 = graph.getFollowing(user2);

        unordered_set<string> set1(following1.begin(), following1.end());
        unordered_set<string> set2(following2.begin(), following2.end());

        int intersection = 0;
        for (const auto& user : set1) {
            if (set2.find(user) != set2.end()) {
                intersection++;
            }
        }

        int unionSize = set1.size() + set2.size() - intersection;
        return unionSize == 0 ? 0.0 : static_cast<double>(intersection) / unionSize;
    }

    // Friend recommendation based on Jaccard similarity
    vector<pair<string, double>> recommendFriends(const string& user, int topN = 5) {
        vector<pair<string, double>> recommendations;
        const auto& userFollowing = graph.getFollowing(user);
        unordered_set<string> userFollowingSet(userFollowing.begin(), userFollowing.end());

        for (const auto& potentialFriend : graph.getAllUsers()) {
            // Skip self and already followed users
            if (potentialFriend != user && userFollowingSet.find(potentialFriend) == userFollowingSet.end()) {
                double similarity = jaccardSimilarity(user, potentialFriend);
                if (similarity > 0) {
                    recommendations.emplace_back(potentialFriend, similarity);
                }
            }
        }

        // Sort by similarity (descending)
        sort(recommendations.begin(), recommendations.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        if (topN > 0 && topN < recommendations.size()) {
            recommendations.resize(topN);
        }

        return recommendations;
    }

    // Find top influencers using PageRank
    vector<pair<string, double>> findTopInfluencers(int topN = 5) {
        auto pageRanks = calculatePageRank();
        vector<pair<string, double>> influencers;
        for (const auto& entry : pageRanks) {
            influencers.emplace_back(entry.first, entry.second);
        }

        // Sort by PageRank (descending)
        sort(influencers.begin(), influencers.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        if (topN > 0 && topN < influencers.size()) {
            influencers.resize(topN);
        }

        return influencers;
    }

    // New method to detect and display communities
   void analyzeCommunities() {
    cout << "Starting community detection..." << endl;
    CommunityDetector detector(graph);
    auto communities = detector.detectCommunities(20); // Limit to 20 iterations
    
    cout << "\n=== Detected Communities ===" << endl;
    for (const auto& comm_entry : communities) {
        int communityId = comm_entry.first;
        const auto& members = comm_entry.second;
        cout << "Community " << communityId << " (" << members.size() << " members):\n  ";
        for (size_t i = 0; i < members.size(); ++i) {
            cout << members[i];
            if (i != members.size() - 1) cout << ", ";
            if (i > 0 && i % 5 == 0) cout << "\n  ";
        }
        cout << endl << endl;
    }
}
};

int main() {
    SocialGraph socialNetwork;

    // Build the social network
    socialNetwork.addConnection("Alice", "Bob");
    socialNetwork.addConnection("Alice", "Charlie");
    socialNetwork.addConnection("Bob", "Charlie");
    socialNetwork.addConnection("Charlie", "Alice");  // Mutual connection
    socialNetwork.addConnection("David", "Alice");
    socialNetwork.addConnection("David", "Bob");
    socialNetwork.addConnection("Eve", "Alice");
    socialNetwork.addConnection("Eve", "David");
    socialNetwork.addConnection("Frank", "Eve");
    socialNetwork.addConnection("Grace", "Frank");
    socialNetwork.addConnection("Grace", "Alice");
    socialNetwork.addConnection("Hank", "Grace");
    socialNetwork.addConnection("Hank", "Ivy");
    socialNetwork.addConnection("Ivy", "Hank");
    socialNetwork.addConnection("Ivy", "Jack");
    socialNetwork.addConnection("Jack", "Ivy");

    SocialAnalyzer analyzer(socialNetwork);

    // Influencer Detection
    cout << "=== Top Influencers (PageRank) ===" << endl;
    auto influencers = analyzer.findTopInfluencers();
    for (const auto& entry : influencers) {
        cout << fixed << setprecision(4);
        cout << entry.first << ": " << entry.second << endl;
    }

    // Friend Recommendation
    string currentUser = "Alice";
    cout << "\n=== Friend Recommendations for " << currentUser << " (Jaccard Similarity) ===" << endl;
    auto recommendations = analyzer.recommendFriends(currentUser);
    for (const auto& entry : recommendations) {
        cout << fixed << setprecision(4);
        cout << entry.first << ": " << entry.second << " similarity" << endl;
    }

    // Community Detection
    analyzer.analyzeCommunities();

    return 0;
}