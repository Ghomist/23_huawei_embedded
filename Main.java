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

class Edge {
    public final int id;
    public final int n1;
    public final int n2;
    public final int dist;
    private final boolean[] pass;
    private int cnt = 0;

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

public class Main {
    static final int EDGE_COST = 1000000;
    static final int PASS_COST = 6;

    static int nodeCnt, edgeCnt, transCnt, passCnt, maxDist;
    static int externCnt = 0;
    static Node[] nodes;
    static Edge[] edges;
    static int[][] trans;
    static int[][] dist;
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
        for (int i = 0; i < nodeCnt; ++i) {
            int start = i;
            int end = i;
            for (Edge edge : nodes[start].edges) {
                end = edge.other(start);
                // var min = nodes[start].getMinEdgeTo(end).dist; // 这个地方接口换了
                var min = Util.shortestEdgeBetween(nodes[start], nodes[end]).dist;
                dist[start][end] = min;
                dist[end][start] = min;
            }
        }

        for (int i = 0; i < nodeCnt; ++i) {
            boolean[] visited = new boolean[nodeCnt];
            int start = i;
            int end = i;
            visited[start] = true;
            int minx = Integer.MAX_VALUE;

            for (int j = 0; j < nodeCnt; ++j) {
                minx = Integer.MAX_VALUE;
                for (int k = 0; k < nodeCnt; ++k) {
                    if (!visited[k] && minx > dist[start][k]) {
                        minx = dist[start][k];
                        end = k;
                    }
                }

                visited[end] = true;

                for (int p = 0; p < nodeCnt; ++p) {
                    if (visited[p])
                        continue;
                    int node2 = 0;
                    int distance = Integer.MAX_VALUE;
                    for (Edge edge : nodes[end].edges) {
                        node2 = edge.other(end);
                        // var min = nodes[end].getMinEdgeTo(node2).dist;
                        var min = Util.shortestEdgeBetween(nodes[end], nodes[node2]).dist;
                        if (node2 == p)
                            distance = min;
                    }
                    if (distance != Integer.MAX_VALUE && dist[start][p] > dist[start][end] + distance) {
                        dist[start][p] = dist[start][end] + distance;
                    }
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

    static void testCapacity() {
        var capacity = new int[edgeCnt];
        for (var t : trans) {
            testPath(t[0], t[1], capacity);
        }
        for (int i = 0; i < edgeCnt; ++i) {
            while (capacity[i] > passCnt) {
                // new edge
                var e = edges[i];
                var newD = Util.shortestEdgeBetween(nodes[e.n1], nodes[e.n2]).dist;
                var newE = new Edge(edgeCnt + externCnt, e.n1, e.n2, newD, passCnt);
                externCnt++;
                externEdges.add(newE);
                nodes[e.n1].edges.add(newE);
                nodes[e.n2].edges.add(newE);
                // reduce capacity
                capacity[i] -= passCnt;
            }
        }
    }

    static void testPath(final int start, final int end, final int[] capacity) {
        var preEdge = new Edge[nodeCnt];
        var cost = new int[nodeCnt];
        var visiting = new PriorityQueue<Node>((a, b) -> {
            // compare F ( F = g + h )
            return Integer.compare(cost[a.id] + dist[a.id][end], cost[b.id] + dist[b.id][end]);
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
                int newCost = cost[node.id] + edge.dist;
                // reach max capacity
                if (capacity[edge.id] >= passCnt)
                    newCost += EDGE_COST * (capacity[edge.id] / passCnt);
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
        if (foundEnd) {
            while (preEdge[n] != null) {
                var e = preEdge[n];
                // arise capacity
                capacity[e.id]++;
                // next node
                n = e.other(n);
            }
        }
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

    /**
     * A* 寻路
     * 
     * @param start 起点
     * @param end   终点
     * @param pass  使用的通道
     */
    static Path findPath(final int start, final int end, final int pass) {
        var preEdge = new Edge[nodeCnt];
        var cost = new int[nodeCnt];
        var visiting = new PriorityQueue<Node>((a, b) -> {
            // compare F ( F = g + h )
            return Integer.compare(cost[a.id] + dist[a.id][end], cost[b.id] + dist[b.id][end]);
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
                int newCost = cost[node.id] + edge.dist;

                // pass has been taken, increase the cost
                if (edge.isTaken(pass)) {
                    if (bridges.contains(edge.id)) {
                        // edge is bridge
                        newCost += EDGE_COST / 20;
                    } else {
                        newCost += EDGE_COST * (passCnt - edge.getTakenPassCnt() + 1); // TODO: 以权重控制加边不是最佳方案
                    }
                }
                // increase cost by pass taken
                // newCost += edge.getTakenPassCnt() * PASS_COST;
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
        // bridges = findBridges();
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
