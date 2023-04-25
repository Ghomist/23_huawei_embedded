import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.Set;

import javax.naming.spi.DirStateFactory.Result;

class Node {
    public final int id;
    public final LinkedList<Edge> edges = new LinkedList<>();

    Node(int id) {
        this.id = id;
    }

    @Override
    public String toString() {
        return "(" + id + ")";
    }
}

class NodeForSort {
    // private final Set<Integer> visited = new HashSet<>();
    public final List<Edge> edges = new LinkedList<>();
    // final Edge e;
    final Node n;
    final int g;
    // final int f;

    NodeForSort(Node n, NodeForSort pre, Edge e, int g) {
        // this.e = e;
        this.n = n;
        this.g = g;
        // this.f = f;
        if (pre != null) {
            // this.visited.addAll(pre.visited);
            this.edges.addAll(pre.edges);
        }
        // this.visited.add(n.id);
        this.edges.add(e);
        // f = g + dist[n.id][end];
    }

    public boolean hasVisited(int nodeID) {
        // return visited.contains(nodeID);
        for (var edge : edges) {
            if (edge != null)
                if (edge.n1 == nodeID || edge.n2 == nodeID)
                    return true;
        }
        return false;
    }
}

class Edge {
    public final int id;
    public final int n1;
    public final int n2;
    public final int dist;
    private final boolean[] pass;

    public long flow = 0;

    Edge(int id, int n1, int n2, int dist, int passCnt) {
        this.dist = dist;
        this.id = id;
        this.n1 = n1;
        this.n2 = n2;
        this.pass = new boolean[passCnt];
    }

    public int other(int node1) {
        return node1 == n1 ? n2 : n1;
    }

    private int cnt = 0;

    public void takePass(int pass) {
        this.pass[pass] = true;
        cnt++;
    }

    public boolean isTaken(int pass) {
        return this.pass[pass];
    }

    public int getTakenPassCnt() {
        return cnt;
    }

    @Override
    public String toString() {
        return "[" + n1 + "]~[" + n2 + "](" + dist + ")";
    }
}

class Path {
    final List<Edge> edges;
    final List<Edge> extern;
    final int start;
    final int end;
    final int cost;
    final int dist;
    final int pass;

    Path(int start, int end, List<Edge> edges, List<Edge> extern, int cost, int dist, int pass) {
        this.extern = extern;
        this.end = end;
        this.start = start;
        this.edges = edges;
        this.cost = cost;
        this.dist = dist;
        this.pass = pass;
    }
}

class Util {
    public static int findMinIdx(int[] arr) {
        int min = 0;
        for (int i = 1; i < arr.length; ++i) {
            if (arr[i] < arr[min])
                min = i;
        }
        return min;
    }

    public static long greatThanZero(long n) {
        return n >= 0 ? n : 0;
    }

    public static Edge shortestEdgeBetween(Node a, Node b) {
        Edge res = null;
        for (Edge edge : a.edges) {
            // is edge between
            if (edge.other(a.id) == b.id) {
                // is shorter edge
                if (res == null || edge.dist < res.dist) {
                    res = edge;
                }
            }
        }
        return res;
    }
}

class Solution {
    final Path[] paths;
    final int cost;
    final Set<Constraint> constraints;

    Solution(Path[] paths) {
        this.paths = paths;
        int cost = 0;
        for (var p : paths) {
            cost += p.cost;
        }
        this.cost = cost;
        this.constraints = new HashSet<>();
    }

    Solution(Path[] paths, Set<Constraint> constraints) {
        this.paths = paths;
        int cost = 0;
        for (var p : paths) {
            cost += p.cost;
        }
        this.cost = cost;
        this.constraints = constraints;
    }

}

class Constraint {
    final int edgeID;
    final int pass;

    public Constraint(int edgeID, int pass) {
        this.edgeID = edgeID;
        this.pass = pass;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || !(obj instanceof Constraint))
            return false;
        var other = (Constraint) obj;
        return edgeID == other.edgeID && pass == other.pass;
    }

    @Override
    public int hashCode() {
        return Objects.hash(edgeID, pass);
    }
}

// class PassCondition {
// Edge edge;// every edge
// boolean[][] pass;// every edge's every pass's condition

// PassCondition(Edge edge, boolean[][] pass) {
// this.edge = edge;
// this.pass = pass;
// }
// }

public class Main {
    static final long EDGE_COST = 1000000;
    static final int PASS_COST = 6;
    static final int TOTAL_FLOW = 64;
    static final float FLOW_DECADE = 2f;

    static final int K = 6;

    // static final Runtime runtime = Runtime.getRuntime();

    static int nodeCnt, edgeCnt, transCnt, passCnt, maxDist;
    static int externCnt = 0;
    static Node[] nodes;
    static Edge[] edges;
    public static int[][] trans;
    static int[][] dist;
    static long maxFlow = 0;
    static Set<Integer> bridges = new HashSet<>();
    static final List<Edge> externEdges = new ArrayList<>();
    static final Map<Integer, String> result = new HashMap<>();

    static void buildMap() {
        final Scanner in = new Scanner(System.in);
        nodeCnt = in.nextInt();
        edgeCnt = in.nextInt();
        transCnt = in.nextInt();
        passCnt = in.nextInt();
        maxDist = in.nextInt();

        nodes = new Node[nodeCnt];
        for (int i = 0; i < nodeCnt; ++i) {
            nodes[i] = new Node(i);
        }

        edges = new Edge[edgeCnt];
        for (int i = 0; i < edgeCnt; ++i) {
            int n1 = in.nextInt();
            int n2 = in.nextInt();
            int d = in.nextInt();
            var e = new Edge(i, n1, n2, d, passCnt);
            edges[i] = e;
            nodes[n1].edges.add(e);
            nodes[n2].edges.add(e);
        }

        trans = new int[transCnt][3];
        for (int i = 0; i < transCnt; ++i) {
            trans[i][0] = in.nextInt();
            trans[i][1] = in.nextInt();
            trans[i][2] = i;
        }

        dist = new int[nodeCnt][nodeCnt];
        for (var line : dist) {
            Arrays.fill(line, Integer.MAX_VALUE);
        }

        in.close();
    }

    /**
     * 生成距离场
     */
    static void buildDistField() {
        // ( >u< )丿✨ hi~~ [❤ love u, my dear ❤] ~~
        // Zzq ❤ Wxr
        for (int i = 0; i < dist.length; ++i) {
            dist[i][i] = 0;
        }
        // for (int i = 0; i < nodeCnt; ++i) {
        // int start = i;
        // int end = i;
        // for (Edge edge : nodes[start].edges) {
        // end = edge.other(start);
        // // var min = nodes[start].getMinEdgeTo(end).dist; // 这个地方接口换了
        // var min = Util.shortestEdgeBetween(nodes[start], nodes[end]).dist;
        // dist[start][end] = min;
        // dist[end][start] = min;
        // }
        // }

        for (int start = 0; start < nodeCnt; ++start) {
            var visited = new boolean[nodeCnt];
            // for (int i = 0; i < start; ++i) {
            // visited[i] = true;
            // dist[start][i] = dist[i][start];
            // }

            visited[start] = true;

            for (int j = 0; j < nodeCnt; ++j) {
                int node = start;
                int minDist = Integer.MAX_VALUE;
                for (int k = 0; k < nodeCnt; ++k) {
                    if (!visited[k] && minDist > dist[start][k]) {
                        minDist = dist[start][k];
                        node = k;
                    }
                }

                visited[node] = true;

                for (var edge : nodes[node].edges) {
                    int next = edge.other(node);

                    if (visited[next])
                        continue;

                    int d = dist[start][node] + edge.dist + 1;
                    if (d < dist[start][next])
                        dist[start][next] = d;
                }

                // for (int p = 0; p < nodeCnt; ++p) {
                // if (visited[p])
                // continue;
                // int node2 = 0;
                // int distance = Integer.MAX_VALUE;
                // for (Edge edge : nodes[node].edges) {
                // node2 = edge.other(node);
                // // var min = nodes[end].getMinEdgeTo(node2).dist;
                // var min = Util.shortestEdgeBetween(nodes[node], nodes[node2]).dist;
                // if (node2 == p)
                // distance = min;
                // }
                // if (distance != Integer.MAX_VALUE && dist[start][p] > dist[start][node] +
                // distance) {
                // dist[start][p] = dist[start][node] + distance;
                // }
                // }
            }
        }
    }

    static void buildFlow() {
        for (var t : trans) {
            buildFlow(t[0], t[1]);
        }
        for (var e : edges) {
            if (e.flow > maxFlow) {
                maxFlow = e.flow;
            }
        }
    }

    static void buildFlow(final int start, final int end) {
        var cnt = new int[nodeCnt];
        // int remainFlow = TOTAL_FLOW;
        int remainK = K;
        var visiting = new PriorityQueue<NodeForSort>(
                (a, b) -> Integer.compare(a.g + dist[a.n.id][end], b.g + dist[b.n.id][end]));
        visiting.add(new NodeForSort(nodes[start], null, null, 0));
        while (!visiting.isEmpty() && remainK > 0) {
            var node = visiting.poll();
            cnt[node.n.id]++;
            // find end
            if (node.n.id == end) {
                for (var e : node.edges) {
                    if (e != null)
                        e.flow += 1;
                }
                // remainFlow /= FLOW_DECADE;
                remainK--;
                continue;
            }
            // ? jias you
            if (cnt[node.n.id] > 7)
                continue;
            // traverse all edges
            for (var edge : node.n.edges) {
                // next node idx
                int next = edge.other(node.n.id);
                // expend new node if current path hasn't been to next node
                if (!node.hasVisited(next)) {
                    visiting.add(new NodeForSort(nodes[next], node, edge, node.g + edge.dist));
                }
            }
        }
    }

    static void sortTrans() {
        var map = new HashMap<Integer, Integer>();
        for (var t : trans) {
            var key = Objects.hash(t[0], t[1]);
            // var key = t[0] + " " + t[1];
            if (map.containsKey(key)) {
                map.put(key, map.get(key) + 1);
            } else {
                map.put(key, 1);// 1204034
            }
        }
        Arrays.sort(trans, (a, b) -> {
            // var keyA = a[0] + " " + a[1];
            // var keyB = b[0] + " " + b[1];
            var keyA = Objects.hash(a[0], a[1]);
            var keyB = Objects.hash(b[0], b[1]);
            int cmp = -Integer.compare(map.get(keyA), map.get(keyB));
            if (cmp != 0)
                return cmp;
            return -Integer.compare(dist[a[0]][a[1]], dist[b[0]][b[1]]);
        });
    }

    static Set<Integer> findBridges() {
        var ansList = new HashSet<Integer>();
        dfs(ansList, new int[nodeCnt], new int[nodeCnt], 0, -1, 1);
        return ansList;
    }

    static void dfs(Set<Integer> ans, int[] dfn, int[] low, int index, int parent, int cnt) {
        dfn[index] = low[index] = cnt;
        for (var edge : nodes[index].edges) {
            int nxt = edge.other(index);
            if (nxt == parent)
                continue;
            if (dfn[nxt] != 0) {
                // visited
                if (dfn[nxt] < low[index]) {
                    low[index] = dfn[nxt];
                }
            } else {
                dfs(ans, dfn, low, nxt, index, cnt + 1);
                if (low[nxt] < low[index]) {
                    low[index] = low[nxt];
                }
            }
        }
        for (var edge : nodes[index].edges) {
            int nxt = edge.other(index);
            if (nxt == parent)
                continue;
            if (low[nxt] > dfn[index]) {
                // find a bridge
                ans.add(edge.id);
            }
        }
    }

    // findpath with Constraint(No Condtraint:-1 -1)
    // edge's passNum is nowalk
    static Path findPath(int start, int end, Set<Constraint> constraints) {

        return null;
    }

    static Solution bestSolution() {
        // all solutions queue
        var minCost = new PriorityQueue<Solution>((a, b) -> {
            return Integer.compare(a.cost, b.cost);
        });

        // first solution (no constraint)
        var paths = new Path[transCnt]; // haolihai!!
        for (int i = 0; i < transCnt; i++) {
            paths[i] = findPath(trans[i][0], trans[i][1], null);
        }
        minCost.add(new Solution(paths));

        // dive into finding best solution
        while (true) {
            var sol = minCost.poll();
            paths = sol.paths;
            boolean judge = false;
            for (int i = 0; i < transCnt - 1; i++) {
                // 加油！！加油！！你打算用 哪一个思路呀 不知道 写写看
                for (int j = i + 1; j < transCnt; j++) {
                    for (int p = 0; p < paths[i].edges.size(); p++) {
                        Edge edge = paths[i].edges.get(p);
                        for (int k = 0; k < paths[j].edges.size(); p++) {
                            if (edge == paths[j].edges.get(k) && paths[i].pass == paths[j].pass) {
                                Constraint newC = new Constraint(p, paths[i].pass);
                                Set<Constraint> cSet = new HashSet<>();
                                cSet.add(newC);
                                cSet.addAll(sol.constraints);
                                Path path = findPath(trans[i][0], trans[i][1], cSet);
                                paths[i] = path;
                                Solution result1 = new Solution(paths, cSet);
                                Set<Constraint> lSet = new HashSet<>();
                                lSet.add(newC);
                                lSet.addAll(sol.constraints);
                                path = findPath(trans[j][0], trans[j][1], lSet);
                                paths[j] = path;
                                Solution result2 = new Solution(paths, lSet);
                                minCost.add(result1);
                                minCost.add(result2);
                                judge = true;
                                break;
                            }
                        }
                        if (judge == true)
                            break;
                    }
                    if (judge == true)
                        break;
                }
                if (judge == true)
                    break;
            }
            if (judge == false)
                return sol;
        }

    }

    /**
     * A* 寻路
     * 
     * @param start 起点
     * @param end   终点
     * @param pass  使用的通道
     */
    static Path findPath(final int start, final int end, final int pass) {
        var preEdge = new Edge[nodeCnt];
        var cost = new double[nodeCnt];
        var visiting = new PriorityQueue<Node>((a, b) -> {
            // compare F ( F = g + h )
            return Double.compare(cost[a.id] + dist[a.id][end], cost[b.id] + dist[b.id][end]);
        });
        Arrays.fill(cost, -1);
        cost[start] = 0;

        boolean foundEnd = false;
        visiting.add(nodes[start]);
        while (!visiting.isEmpty()) {
            var node = visiting.poll();
            if (node.id == end) {
                foundEnd = true;
                break;
            }
            // traverse all edges
            for (var edge : node.edges) {
                // next node idx
                int next = edge.other(node.id);
                // cal new g
                double newCost = cost[node.id] + (double) edge.dist / maxDist + 1;
                // pass has been taken, increase the cost
                if (edge.isTaken(pass)) {
                    if (bridges.contains(edge.id)) {
                        // edge is bridge
                        newCost += EDGE_COST * (2 - ((double) edge.flow / maxFlow)) / 20;
                    } else {
                        newCost += EDGE_COST * (2 - ((double) edge.flow / maxFlow));
                        // * (passCnt + edge.getTakenPassCnt() + 1);
                    }
                }
                // update g
                if (cost[next] == -1 || newCost < cost[next]) {
                    // update array
                    cost[next] = newCost;
                    // update pre node
                    preEdge[next] = edge;
                    // set new node to visiting (do NOT repeat old node!)
                    if (!visiting.contains(nodes[next]))
                        visiting.add(nodes[next]);
                }
            }
        }
        int n = end;
        int externId = externCnt + edgeCnt;
        if (foundEnd) {
            var extern = new LinkedList<Edge>();
            var path = new LinkedList<Edge>();
            int dist = 0;
            while (preEdge[n] != null) {
                var e = preEdge[n];
                // if need new edge
                if (e.isTaken(pass)) {
                    var newD = Util.shortestEdgeBetween(nodes[e.n1], nodes[e.n2]).dist;
                    var newE = new Edge(externId, e.n1, e.n2, newD, passCnt);
                    externId++;
                    extern.add(newE);
                    e = newE;
                }
                // make path
                path.addFirst(e);
                // cal total dist
                dist += e.dist;
                // next node
                n = e.other(n);
            }
            return new Path(start, end, path, extern, cost[end], dist, pass);
        } else {
            return null;
        }
    }

    static List<Integer> setAmp(Path path) {
        var ans = new LinkedList<Integer>();
        int crt = path.start;
        int edgeLength = 0;
        for (Edge e : path.edges) {
            int nxt = e.other(crt);
            edgeLength += e.dist;
            if (edgeLength > maxDist) {
                ans.add(crt);
                edgeLength = e.dist;
            } else if (edgeLength == maxDist) {
                if (nxt != path.end) // 防止放大器不在路径上
                    ans.add(nxt);
                edgeLength = 0;
            }
            crt = nxt;
        }
        return ans;
    }

    static void applyPath(Path path, int index) {
        var edges = path.edges;
        var pass = path.pass;

        var ampList = setAmp(path);

        var sb = new StringBuilder();
        sb.append(pass).append(' ')
                .append(edges.size()).append(' ')
                .append(ampList.size()).append(' ');

        // apply extern edges
        externEdges.addAll(path.extern);
        externCnt += path.extern.size();
        for (var ex : path.extern) {
            nodes[ex.n1].edges.add(ex);
            nodes[ex.n2].edges.add(ex);
        }

        // edge
        for (var edge : edges) {
            // take pass
            edge.takePass(pass);
            // output
            sb.append(edge.id).append(' ');
        }

        // make amplifier
        for (var amp : ampList) {
            sb.append(amp).append(' ');
        }

        sb.deleteCharAt(sb.length() - 1);
        // System.out.println(sb);
        result.put(index, sb.toString());
    }

    static void output() {
        System.out.println(externCnt);
        for (var ex : externEdges) {
            System.out.println(ex.n1 + " " + ex.n2);
        }
        // plz!!!!
        for (int i = 0; i < transCnt; i++) {
            var resultI = result.get(i);
            System.out.println(resultI);
        }
        // for (var line : result) {
        // System.out.println(line);
        // }
    }

    public static void main(String[] args) {
        buildMap();
        buildDistField();
        // buildFlow();
        bridges = findBridges();
        // sortTrans();
        // testCapacity();
        for (var t : trans) {
            Path minPath = null;
            for (int pass = 0; pass < passCnt; pass++) {
                var path = findPath(t[0], t[1], pass);
                if (path != null) {
                    int cost = path.cost;
                    if (minPath == null || cost < minPath.cost) {
                        minPath = path;
                    }
                }
            }
            applyPath(minPath, t[2]);
        }
        output();
    }

}
