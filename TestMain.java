public class TestMain {
    public static void main(String[] args) {
        int n = 64;
        int cnt = 0;
        double d = 2;
        while (n > 0) {
            n /= d;
            cnt++;
        }
        System.out.println(cnt);
    }
}
