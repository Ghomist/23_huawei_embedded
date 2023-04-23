import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

public class TestGenerator {

    private static final int PASS_CNT = 80;
    public static final int NODE_CNT = 1000 ;
    public static final int TRANS_CNT = 8000;

    public static final int MAX_DIST = 100;
    public static final int MAX_EDGE_CNT = 4000;
    public static final double EX_EDGE_RATE = 0.3; // 随机添加新边的概率
    public static final double DUP_EDGE_RATE = 0.1; // 随机添加多重边的概率
    public static final double MAX_DUP_CNT = 3; // 多重边的最大数量

    public static final Random RAND = new Random();

    public static int[][] generateGraph(int n) {
        // 初始化邻接矩阵
        int[][] adjMatrix = new int[n][n];

        // 初始化边的列表
        System.err.println("Init edges");
        List<int[]> edges = new LinkedList<>();
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                edges.add(new int[] { i, j, RAND.nextInt(n * n) });
                // System.err.print("\rInit edges: " + i + ',' + j);
            }
        }

        // 按照权值从小到大排序边的列表
        // Collections.shuffle(edges);
        edges.sort(Comparator.comparingInt(e -> e[2]));

        System.err.println("Gen edges");
        // Kruskal 算法生成连通图
        int edgeCnt = 0;
        int[] parent = new int[n];
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }
        for (int e = 0; e < edges.size(); ++e) {
            var edge = edges.get(e);
            int u = find(parent, edge[0]);
            int v = find(parent, edge[1]);
            if (u != v) {
                adjMatrix[edge[0]][edge[1]] = edge[2];
                adjMatrix[edge[1]][edge[0]] = edge[2];
                parent[u] = v;
                edgeCnt++;
            }
            // System.err.print("\rGen edges: " + e + " of " + edgeCnt);
        }

        System.err.println();

        // 随机添加额外的边
        System.err.println("Gen ex edge");
        Collections.shuffle(edges);
        for (int e = 0; e < edges.size(); ++e) {
            var edge = edges.get(e);
            if (edgeCnt > MAX_EDGE_CNT) {
                break;
            }
            int u = edge[0];
            int v = edge[1];
            if (adjMatrix[u][v] == 0 && RAND.nextDouble() < EX_EDGE_RATE) {
                adjMatrix[u][v] = edge[2];
                adjMatrix[v][u] = edge[2];
                edgeCnt++;
            }
            // System.err.print("\rGen ex: " + e + " of " + edgeCnt);
        }

        return adjMatrix;
    }

    // 并查集查找节点所在的连通分量
    private static int find(int[] parent, int node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        return node;
    }

    public static List<int[]> adjMatrixToList(int[][] adjMatrix) {
        System.err.println("To list");
        List<int[]> edgeList = new ArrayList<>();
        int n = adjMatrix.length;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (adjMatrix[i][j] > 0) {
                    // 随机化所有的边长
                    edgeList.add(new int[] { i, j, RAND.nextInt(MAX_DIST) });
                    // 随机添加多重边
                    for (int dupI = 0; dupI < MAX_DUP_CNT; dupI++) {
                        if (RAND.nextDouble() < DUP_EDGE_RATE) {
                            edgeList.add(new int[] { i, j, RAND.nextInt(MAX_DIST) });
                        }
                    }
                }
            }
        }
        return edgeList;
    }

    public static void main(String[] args) throws FileNotFoundException {
        PrintStream ps = new PrintStream(new FileOutputStream("sample.txt"));
        System.setOut(ps);

        List<int[]> edgeList = adjMatrixToList(generateGraph(NODE_CNT));
        System.out.println(NODE_CNT + " " + edgeList.size() + " " + TRANS_CNT + " " + PASS_CNT + " " + MAX_DIST);
        for (int[] edge : edgeList) {
            System.out.println(edge[0] + " " + edge[1] + " " + edge[2]);
        }
        System.err.println("Gen trans");
        for (int t = 0; t < TRANS_CNT; t++) {
            int from = RAND.nextInt(NODE_CNT);
            int to = RAND.nextInt(NODE_CNT);
            while (from == to) {
                to = RAND.nextInt(NODE_CNT);
            }
            System.out.println(from + " " + to);
        }

        ps.close();
    }
}
