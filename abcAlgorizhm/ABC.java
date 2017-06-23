package abcAlgorizhm;

public class ABC {

	private int honeyNum;
	private int tryLimit;
	private int step;
	private int allStep;
	private int outputRange;
	private boolean stop;

	private Honey[] honeys;
	static int num;
	static int n;

	public ABC(int honeyNum, int tryLimit, int outputRange, int num, int n) {
		this.step = 0;
		this.allStep = 0;
		this.stop = false;
		this.honeyNum = honeyNum;
		this.tryLimit = tryLimit;
		this.outputRange = outputRange;
		ABC.num = num;
		ABC.n = n;
		getInitHoneys();
	}

	public void work() {
		while (!stop) {
			step++;
			employedBeeFase();
			onlookerBeeFase();
			if (step == outputRange) {
				step = 0;
				allStep += outputRange;
				System.out.println(allStep + "\t" + Honey.min);
			}
			if (allStep == outputRange * 100) {
				System.out.println("最適値:" + Honey.min);
				System.out.print("x=(");
				for (int j = 0; j < n; j++) {
					if (j != n - 1) {
						System.out.print(Honey.minX[j] + ",");
					} else {
						System.out.print(Honey.minX[j] + ")\n");
					}
				}
				stop = true;
			}
			scoutBeeFase();
		}
		printResult();
	}

	private void employedBeeFase() {
		for (int i = 0; i < honeyNum; i++) {
			int anotherHoney = searchAnotherHoney(i);
			moveToAnother(honeys[i], honeys[anotherHoney]);
		}
	}

	private void onlookerBeeFase() {
		for (int i = 0; i < honeyNum; i++) {
			int change = 0;
			int A = (int) (honeyNum * Math.random());
			int B = (int) (honeyNum * Math.random());
			int C = (int) (honeyNum * Math.random());
			int D = (int) (honeyNum * Math.random());
			int E = (int) (honeyNum * Math.random());
			if (honeys[A].fitness > honeys[B].fitness) {
				change = A;
			} else {
				change = B;
			}
			if (honeys[change].fitness < honeys[C].fitness) {
				change = C;
			}
			if (honeys[change].fitness < honeys[D].fitness) {
				change = D;
			}
			if (honeys[change].fitness < honeys[E].fitness) {
				change = E;
			}
			int anotherHoney = searchAnotherHoney(change);
			honeys[i]=moveToAnother(honeys[change], honeys[anotherHoney]);
		}
	}

	private void scoutBeeFase() {
		for (int i = 0; i < honeyNum; i++) {
			if (honeys[i].tryCounter >= this.tryLimit) {
				honeys[i] = new Honey(num, n);
			}
		}
	}

	private void getInitHoneys() {
		honeys = new Honey[honeyNum];
		for (int i = 0; i < honeyNum; i++) {
			honeys[i] = new Honey(num, n);
		}
	}

	private int searchAnotherHoney(int i) {
		boolean check = false;
		int anotherHoney = 0;
		while (!check) {
			anotherHoney = (int) (honeyNum * Math.random());
			if (anotherHoney == i || anotherHoney == honeyNum)
				check = false;
			else
				check = true;
		}
		return anotherHoney;
	}

	private Honey moveToAnother(Honey baseHoney, Honey toHoney) {
		Honey newHoney = new Honey(n);
		for (int j = 0; j < n; j++) {
			boolean check = false;
			double x = 0;
			while (!check) {
				double r = 2 * (Math.random() - 0.5);
				x = baseHoney.x[j] + r * (baseHoney.x[j] - toHoney.x[j]);
				check = Honey.isRange(x);
			}
			newHoney.x[j] = x;
		}
		newHoney.calcurateFitness();
		if (baseHoney.fitness < newHoney.fitness)
			return newHoney;
		else{
			baseHoney.tryCounter++;
			return baseHoney;
		}
	}

	private void printResult() {
		// TODO 自動生成されたメソッド・スタブ

	}
}
