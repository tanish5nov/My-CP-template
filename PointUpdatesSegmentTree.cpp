class SegmentTree {
public:

	/*
		al -> array left index for querying
		ar -> array right index for querying
		pl -> provided left index of range to query
		pr -> provided right index of range to query
		index -> current index in tree
		updIndex -> index to be updated in point updates
		updVal -> new value to update that index with in point updates
	*/

	vector<int>tree;
	void init(int sz) {
		int twoPower = 1;
		while (twoPower < sz) {
			twoPower *= 2;
		}
		tree.resize(twoPower * 2ll);
	}

	int buildTree(int al, int ar, int index, vector<int>&arr) {
		if (al > ar) {
			return 0;
		}
		if (al == ar) {
			tree[index] = arr[al];
			return arr[al];
		}
		int mid = (al + ar) / 2;
		int left = buildTree(al, mid, (2 * index) + 1, arr);
		int right = buildTree(mid + 1, ar, (2 * index) + 2, arr);
		return tree[index] = left + right;
	}

	int queryPoint(int pl, int pr, int al, int ar, int index) {
		if (al > ar || ar < pl || al > pr) {
			return 0;
		}
		if (pl <= al && ar <= pr) {
			// bug(pl, pr, al, ar);
			return tree[index];
		}
		int mid = (al + ar) / 2;
		int left = queryPoint(pl, pr, al, mid, (2 * index) + 1);
		int right = queryPoint(pl, pr, mid + 1, ar, (2 * index) + 2);
		return left + right;
	}

	void pointUpdates(int al, int ar, int updIndex, int updVal, int index, vector<int>&arr) {
		if (al > ar) {
			return;
		}

		if (al == ar && al == updIndex) {
			arr[updIndex] = updVal;
			tree[index] = updVal;
			return;
		}

		int mid = (al + ar) / 2;
		if (updIndex <= mid) {
			pointUpdates(al, mid, updIndex, updVal, (2 * index) + 1, arr);
		} else {
			pointUpdates(mid + 1, ar, updIndex, updVal, (2 * index) + 2, arr);
		}
		int left = tree[(2 * index) + 1];
		int right = tree[(2 * index) + 2];
		tree[index] = left + right;
	}
};
