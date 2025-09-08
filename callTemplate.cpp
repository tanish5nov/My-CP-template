#include<bits/stdc++.h>
#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
#ifndef ONLINE_JUDGE
#pragma comment(linker, "/STACK:1073741824")
#endif


//------------------------------ TEMPLATE STARTS -------------------------//
#define int long long
template<class T> using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<class T> using ordered_multiset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;
#ifndef ONLINE_JUDGE
#define bug(...) __f (#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>void __f (const char* name, Arg1&& arg1) { cout << name << " : " << arg1 << endl; }
template <typename Arg1, typename... Args>void __f (const char* names, Arg1&& arg1, Args&&... args) {const char* comma = strchr (names + 1, ','); cout.write (names, comma - names) << " : " << arg1 << " | "; __f (comma + 1, args...);}
#define pr(a) {cout<<#a<<": ";for(auto it: a)cout<<it<<" ";cout<<endl;}
#define pm(a) {cout<<#a<<": "<<endl;for(auto it: a)cout<<it.first<<" "<<it.second<<" "<<endl;}
#else
#define bug(...) __f (#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>void __f (const char* name, Arg1&& arg1) {}
template <typename Arg1, typename... Args>void __f (const char* names, Arg1&& arg1, Args&&... args) {}
#define pr(a) {for(auto it: a)cout<<it<<" ";cout<<endl;}
#define pm(a) {for(auto it: a)cout<<it.first<<" "<<it.second<<" "<<endl;}
#endif
#define all(a) a.begin(),a.end()
#define in(a) {for(auto &i: a)cin>>i;}
#define out(a) {for(auto &i: a)cout<<i<<" ";cout<<endl;}
#define out2d(a){for(auto i: a){out(i);}}
#define sort(a) sort(a.begin(),a.end());
#define rev(a) reverse(a.begin(),a.end());
#define endl "\n"
#define pb push_back
#define pp pop_back
#define yes cout<<"YES"<<endl;
#define no cout<<"NO"<<endl;
#define syes cout<<"Yes"<<endl;
#define sno cout<<"No"<<endl;
#define BLUE_AND_PURPLE {ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);}
#define ceil(x,y) (x+y-1)/y



/*--------------------------FACTORIAL - MOD----------------------*/
const int MOD = 1e9 + 7;
int mod(int x, int m) {
    return (x % m + m) % m;
}
int modMul(int a, int b, int m) {
    a = mod(a, m);
    b = mod(b, m);
    return ((a * b) % m);
}
int modAdd(int a, int b, int m) {
    a = mod(a, m);
    b = mod(b, m);
    return ((a + b) % m);
}
int modSub(int a, int b, int m) {
    a = mod(a, m);
    b = mod(b, m);
    return ((a - b + m) % m);
}
int modPow(int x, int y, int m) {
    if (y == 0) {
        return 1;
    }
    int ans1 = modPow(x, y / 2, m);
    int ans2 = ans1;
    if (y & 1) {
        return modMul(ans1, modMul(ans2, x, m), m);
    } else {
        return modMul(ans1, ans2, m);
    }
}
int invMod(int a, int m) {
    //(1/a)%m = a^(m-2)
    return modPow(a, m - 2, m);
}
int modDiv(int a, int b, int m) {
    return modMul(a, invMod(b, m), m);
}

int recursiveNcR(int n, int r, vector<vector<int>>&dp) {
    if (n < r) {
        return dp[n][r] = 0ll;
    }
    if (r == 0) {
        return dp[n][r] = 1ll;
    }
    if (r == 1) {
        return dp[n][r] = n;
    }
    if (n == 1) {
        return dp[n][r] = 1ll;
    }
    if (dp[n][r] != -1)return dp[n][r];
    return dp[n][r] = (recursiveNcR(n - 1, r - 1, dp) + recursiveNcR(n - 1, r, dp)) % MOD;
}

int dereangement(int n, vector<int>&dp) {
    if (n == 0) {
        return dp[n] = 1;
    } else if (n == 1) {
        return dp[n] = 0;
    } else if (n == 2) {
        return dp[n] = 1;
    }
    if (dp[n] != -1)return dp[n];
    int ff = modSub(n, 1ll, MOD);
    int ss = dereangement(n - 1, dp);
    int tt = dereangement(n - 2, dp);
    int sum = modAdd(ss, tt, MOD);
    int mul = modMul(ff, sum, MOD);
    return dp[n] = mul;
}

const int factN = 2e5 + 10;
int fact[factN];
void preFact(int m) {
    fact[0] = 1;
    for (int i = 1; i < factN; i++) {
        fact[i] = modMul(fact[i - 1], i, m);
    }
}

int ncr(int n, int r, int m) {
    if (n < r)return 0;
    int den = modMul(fact[r], fact[n - r], MOD);
    int num = fact[n];
    int ans = modDiv(num, den, m);
    return ans;
}



/*-------------------------- SPARSE- TABLE ----------------------*/
class sparseTable {
private:

    int N, K;
    vector<int>pre;
    vector<vector<int>>dp;
    int f(int ff, int ss) {
        //ONLY CHANGES ARE REQUIRED HERE
        return (ff & ss);
    }

public:

    void init(int n, int k) {
        N = n;
        K = k;
        pre.resize(N + 1);
        dp.resize(N + 1);
        for (int i = 0; i < dp.size(); i++) {
            dp[i].resize(K);
        }
        for (int i = 1; i < N; i++) {
            pre[i] = log2l(i);
        }
    }

    void build(vector<int>arr) {
        for (int i = 0; i < arr.size(); i++) {
            dp[i][0] = arr[i];
        }
        for (int j = 1; j < K; j++) {
            for (int i = 0; i + (1 << j) < N; i++) {
                dp[i][j] = f(dp[i][j - 1], dp[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    int query(int L, int R) {
        // 0 Based Indexing Arguments to be passed
        int j = pre[R - L + 1];
        int ans = f(dp[L][j], dp[R - (1 << j) + 1][j]);
        return ans;
    }
};



/*-------------------------- 2D - SPARSE- TABLE ----------------------*/

class SparseTable2D {
public:
    //0 based indexing
    vector<vector<vector<vector<int>>>>sparse;
    vector<int>log;
    int R, C;

    int f(int X, int Y) {
        return max(X, Y);
    }

    int query(int x1, int y1, int x2, int y2, bool build) {
        int x_sz = x2 - x1 + 1;
        int y_sz = y2 - y1 + 1;
        int k1 = (x_sz == 1) ? 0 : log[build ? (x_sz - 1) : x_sz];
        int k2 = (y_sz == 1) ? 0 : log[build ? (y_sz - 1) : y_sz];
        int NW = sparse[k1][k2][x1][y1];
        int NE = sparse[k1][k2][x1][y_sz - (1 << k2) + y1];
        int SW = sparse[k1][k2][x_sz - (1 << k1) + x1][y1];
        int SE = sparse[k1][k2][x_sz - (1 << k1) + x1][y_sz - (1 << k2) + y1];
        return f(f(NW, NE), f(SW, SE));
    }

    void init(vector<vector<int>>arr, int RR, int CC) {
        R = RR;
        C = CC;
        log.resize(max(R, C) + 1000);

        for (int i = 0; i < log.size(); i++) {
            log[i] = log2l(i);
        }

        int k1 = log[R] + 1;
        int k2 = log[C] + 1;
        sparse.resize(k1);
        for (int i = 0; i < k1; i++) {
            sparse[i].resize(k2);
            for (int j = 0; j < k2; j++) {
                sparse[i][j].resize(R);
                for (int x = 0; x < R; x++) {
                    sparse[i][j][x].resize(C);
                }
            }
        }
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < C; j++) {
                sparse[0][0][i][j] = arr[i][j];
            }
        }

        for (int h = 0; h < k1; h++) {
            for (int v = 0; v < k2; v++) {
                if (!(h == 0 && v == 0)) {
                    for (int i = 0; i + (1 << h) <= R; i++) {
                        for (int j = 0; j + (1 << v) <= C; j++) {
                            sparse[h][v][i][j] = query(i, j, i + (1 << h) - 1, j + (1 << v) - 1, true);
                        }
                    }
                }
            }
        }
    }
};



/*-------------------------- Disjoint Sets ----------------------*/

class DSU {
public:
    vector<int>parent, sizes;
    void init(int n) {
        parent.resize(n + 1);
        sizes.resize(n + 1);
        for (int i = 1; i <= n; i++) {
            parent[i] = i;
            sizes[i] = 1;
        }
    }
    int finder(int node) {
        if (parent[node] == node) {
            return node;
        }
        return parent[node] = finder(parent[node]);
    }

    void merge(int node1, int node2) {
        node1 = finder(node1);
        node2 = finder(node2);

        if (sizes[node1] > sizes[node2]) {
            parent[node2] = node1;
            sizes[node1] += sizes[node2];
        } else {
            parent[node1] = node2;
            sizes[node2] += sizes[node1];
        }
    }
};



/*-------------------------- GCD - LCM ----------------------*/
int gcd(int x, int y) {
    if (y == 0) return x;
    return gcd(y, x % y);
}
int lcm(int x, int y) {
    return ((x * y) / gcd(x, y));
}



/*-------------------------- PRIME SIEVE ----------------------*/
const int maxn = 1e6 + 8;
int spf[maxn];
void sieve() {
    spf[1] = 0;
    spf[0] = 1;
    for (int i = 2; i < maxn; i++) {
        spf[i] = i;
    }
    for (int i = 2; i * i < maxn; i++) {
        if (spf[i] == i) {
            for (int j = i * i; j < maxn; j += i) {
                if (spf[j] == j) {
                    spf[j] = i;
                }
            }
        }
    }
}

int phi[maxn];
void totient() {
    // phi[a * b] = phi[a] * phi[b] given that a,b are co-prime i.e gcd(a,b) == 1
    // sieve is required for this function
    // phi[i] -> number of co-primes of i in range [1, maxn)
    // totient formula f(x) = x * (1 - 1/p1) * (1 - 1/p2) ... *(1 - 1/pn)
    // where p1 to pn represent the prime factorization of x
    for (int i = 1; i < maxn; i++) {
        phi[i] = i;
    }
    for (int i = 2; i < maxn; i++) {
        if (spf[i] == i) {
            for (int j = i; j < maxn; j += i) {
                /*
                    (just for the sake of reminding if I am stuck)
                    explanation of this relation:
                    we first set all numbers phi function to i
                    now we know the formula for totient of a number x if
                    f(x) = x * (1 - 1/p1) * (1 - 1/p2) ... *(1 - 1/pn)
                    which is eventually derived from f(x) = x - (x/p1) - (x/p2) + (x/p1p2)....
                    so we are doing the same thing that for every number which is divided by
                    the prime i do it as x -= (x / pi)
                */
                phi[j] -= phi[j] / i;
            }
        }
    }
}

void getFactors(int n, vector<int>&fcs) {
    while (n != 1) {
        fcs.push_back(spf[n]);
        n /= spf[n];
    }
}

void getDivisors(int n, set<int>&divs) {
    for (int i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            if (n / i == i) {
                divs.insert(i);
            } else {
                divs.insert(n / i);
                divs.insert(i);
            }
        }
    }
}

//------------------------------ TEMPLATE ENDS -------------------------//











void solve(int TestCase) {

}

int32_t main() {
    ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
    // preFact(MOD);
    // sieve();
    // totient();
    int tests = 1;
    cin >> tests;
    for (int i = 1; i <= tests; i++)solve(i);
}











e
