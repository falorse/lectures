package abcAlgorizhm;

import java.util.Date;

public class Tester {

	public void test(int num, int n, int honeyNum, int tryLimit, int outputRange) {

		System.out.println("問題番号:" + num + "\n" + "次元:" + n + "\n" + "探索点の数:"
				+ honeyNum + "\n" + "試行カウンタの限界値:" + tryLimit + "\n" + "出力間隔:"
				+ outputRange);
		Date start = new Date();

		Honey.min = 2560000;
		ABC abc = new ABC(honeyNum, tryLimit, outputRange, num, n);
		abc.work();

		Date finish = new Date();
		System.out.println("計算時間:" + (finish.getTime() - start.getTime())
				+ "(msec)");

	}

	public static void main(String args[]) {
		Tester test = new Tester();
		// 問題番号、次元、探索点の数、
		//試行カウンタの限界値、出力間隔 を入力してください。
		test.test(4, 100, 1000, 100, 100);
	}

}
