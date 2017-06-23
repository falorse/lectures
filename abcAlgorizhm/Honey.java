package abcAlgorizhm;

public class Honey {

	int tryCounter;
	double fitness;
	public double x[];
	static int n;
	static double min;
	static double minX[];
	static int num;

	public Honey(int num,int n){
		this.tryCounter=0;
		this.fitness=0;
		Honey.num=num;
		Honey.n=n;
		this.x=new double[n];
		initHoney();
	}

	public Honey(int n) {
		this.tryCounter=0;
		this.fitness=0;
		Honey.num=num;
		Honey.n=n;
		this.x=new double[n];
		}

	public void initHoney(){
		this.x=setX();
		calcurateFitness();
		this.tryCounter=0;
	}

	void calcurateFitness() {
		double sum = 0;
		int N = x.length;
		if (num == 1) {
			for (int i = 0; i < N; i++) {
				sum += pow2(x[i]);
			}
		} else if (num == 2) {
			sum += 10 * N;
			for (int i = 0; i < N; i++) {
				sum += pow2(x[i]) - 10 * Math.cos(2 * Math.PI * x[i]);
			}
		} else if (num == 3) {
			for (int i = 0; i < N - 1; i++) {
				sum += 100 * pow2(x[i + 1] - pow2(x[i])) + pow2(1 - x[i]);
			}
		} else if (num == 4) {
			sum += 1;
			for (int i = 0; i < N; i++) {
				sum += pow2(x[i]) /4000-  Math.cos(x[i] / Math.sqrt(i+1)) ;
			}
		} else if (num == 5) {
			for (int i = 0; i < N; i++) {
				sum += Math.abs(x[i] * Math.sin(x[i]) + 0.1 * x[i]);
			}
		} else if (num == 6) {
			for (int i = 0; i < N; i++) {
				sum += Math.pow(x[i], 4) - 16 * pow2(x[i]) + 5 * x[i];
			}
		}
		this.fitness=1/(100000+sum);
		if(sum<min){
			Honey.min=sum;
			Honey.minX=this.x;
		}
	}

	private double[] setX() {
		double[] point = new double[n];
		if (num == 1 || num == 2) {
			for (int j = 0; j < n; j++) {
				point[j] = 10 * (Math.random() - 0.5);
			}
		} else if (num == 3) {
			int j = 0;
			while (j < n) {
				double r = 20 * (Math.random() - 0.5);
				if (isRange(r)) {
					point[j] = r;
					j++;
				}
			}
		} else if (num == 4) {
			for (int j = 0; j < n; j++) {
				point[j] = 1200 * (Math.random() - 0.5);
			}
		} else if (num == 5) {
			int j = 0;
			while (j < n) {
				double r = 20 * (Math.random() - 0.5);
				if (isRange(r)) {
					point[j] = r;
					j++;
				}
			}
		} else if (num == 6) {
			int j = 0;
			while (j < n) {
				double r = 10 * (Math.random() - 0.5);
				if (isRange(r)) {
					point[j] = r;
					j++;
				}
			}
		}
		return point;
	}

	public static boolean isRange(double x) {
		if (num == 1 || num == 2) {
			if (-5.0 <= x && x <= 5.0) {
				return true;
			}
		} else if (num == 3) {
			if (-5.0 <= x && x < 10.0) {
				return true;
			}
		} else if (num == 4) {
			if (-600.0 <= x && x <= 600.0) {
				return true;
			}
		} else if (num == 5) {
			if (-10.0 <= x && x < 10.0) {
				return true;
			}
		} else if (num == 6) {
			if (-5.0 <= x && x < 5.0) {
				return true;
			}
		} else {
			return false;
		}
		return false;
	}

	public static double pow2(double x) {
		return x * x;
	}

}
