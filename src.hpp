#ifndef SRC_HPP
#define SRC_HPP

// Intentionally keep includes minimal to match the template expectations.
// The judge environment provides standard headers via its driver and fraction.hpp.
#include "fraction.hpp"

class matrix {
private:
    int m, n;
    fraction **data;

    void allocate(int rows, int cols) {
        m = rows;
        n = cols;
        if (m <= 0 || n <= 0) {
            data = nullptr;
            m = n = 0;
            return;
        }
        data = new fraction *[m];
        for (int i = 0; i < m; ++i) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; ++j) data[i][j] = fraction(0);
        }
    }

public:
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) { allocate(m_, n_); }

    // 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) {
        allocate(obj.m, obj.n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data[i][j] = obj.data[i][j];
            }
        }
    }

    // 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept {
        m = obj.m; n = obj.n; data = obj.data;
        obj.m = obj.n = 0; obj.data = nullptr;
    }

    // 析构函数。
    ~matrix() {
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        data = nullptr; m = n = 0;
    }

    // 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        // Free existing
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        allocate(obj.m, obj.n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) data[i][j] = obj.data[i][j];
        }
        return *this;
    }

    // 重载括号：第i行(1-based)、第j列(0-based)
    fraction &operator()(int i, int j) {
        // i: 1..m, j: 0..n-1
        if (!(1 <= i && i <= m && 0 <= j && j < n)) {
            throw matrix_error();
        }
        return data[i - 1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n == 0 || rhs.m == 0 || lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix res(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; ++i) {
            for (int k = 0; k < lhs.n; ++k) {
                const fraction &lik = lhs.data[i][k];
                if (lik == fraction(0)) continue;
                for (int j = 0; j < rhs.n; ++j) {
                    res.data[i][j] = res.data[i][j] + lik * rhs.data[k][j];
                }
            }
        }
        return res;
    }

    // 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
    matrix transposition() {
        if (m == 0 || n == 0 || data == nullptr) throw matrix_error();
        matrix t(n, m);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                t.data[j][i] = data[i][j];
            }
        }
        return t;
    }

    // 返回矩阵的行列式（高斯消元）。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (m == 0 || n == 0 || m != n || data == nullptr) throw matrix_error();
        // Make a copy
        matrix tmp(*this);
        fraction det(1);
        int sign = 1;
        for (int col = 0; col < n; ++col) {
            int pivot = col;
            // find non-zero pivot
            while (pivot < m && tmp.data[pivot][col] == fraction(0)) ++pivot;
            if (pivot == m) {
                return fraction(0);
            }
            if (pivot != col) {
                // swap rows
                fraction *rowA = tmp.data[col];
                tmp.data[col] = tmp.data[pivot];
                tmp.data[pivot] = rowA;
                sign = -sign;
            }
            fraction piv = tmp.data[col][col];
            // eliminate below
            for (int r = col + 1; r < m; ++r) {
                if (tmp.data[r][col] == fraction(0)) continue;
                fraction factor = tmp.data[r][col] / piv;
                // row_r = row_r - factor * row_col
                for (int c = col; c < n; ++c) {
                    tmp.data[r][c] = tmp.data[r][c] - factor * tmp.data[col][c];
                }
            }
            det = det * piv;
        }
        if (sign == -1) det = fraction(0) - det; // negate
        return det;
    }
};

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction;

    // edges stored implicitly via adjacency + conduction
    // Helper to build Laplacian reduced matrix for nodes 1..n-1
    void build_reduced_conductance(fraction **G) {
        int n = interface_size;
        int k = (n >= 1 ? n - 1 : 0);
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) G[i][j] = fraction(0);
        }
        // For each edge e, endpoints u,v and conductance g
        for (int e = 0; e < connection_size; ++e) {
            // adjacency: row e has +1 at from-1, -1 at to-1
            int u = -1, v = -1;
            for (int j = 0; j < interface_size; ++j) {
                fraction val = adjacency(1 + e, j); // i is 1-based
                if (val == fraction(1)) u = j; // j is 0-based
                else if (val == fraction(0) - fraction(1)) v = j;
            }
            // conductance on diagonal
            fraction g = conduction(1 + e, e); // diagonal entry
            // add contributions
            int iu = u, iv = v;
            if (iu != interface_size - 1) {
                G[iu][iu] = G[iu][iu] + g;
            }
            if (iv != interface_size - 1) {
                G[iv][iv] = G[iv][iv] + g;
            }
            if (iu != interface_size - 1 && iv != interface_size - 1) {
                G[iu][iv] = G[iu][iv] - g;
                G[iv][iu] = G[iv][iu] - g;
            }
        }
    }

    static void solve_linear_system(int n, fraction **A, fraction *b, fraction *x) {
        // Gaussian elimination without pivoting beyond zero checks
        // Forward elimination
        for (int col = 0, row = 0; col < n && row < n; ++col, ++row) {
            int pivot = row;
            while (pivot < n && A[pivot][col] == fraction(0)) ++pivot;
            if (pivot == n) {
                // Singular or no pivot in this column; skip (system is guaranteed solvable)
                --row; // so outer ++row keeps same row index for next column
                continue;
            }
            if (pivot != row) {
                // swap rows in A and b
                fraction *tmp = A[pivot];
                A[pivot] = A[row];
                A[row] = tmp;
                fraction tb = b[pivot]; b[pivot] = b[row]; b[row] = tb;
            }
            fraction piv = A[row][col];
            // eliminate rows below
            for (int r = row + 1; r < n; ++r) {
                if (A[r][col] == fraction(0)) continue;
                fraction factor = A[r][col] / piv;
                for (int c = col; c < n; ++c) {
                    A[r][c] = A[r][c] - factor * A[row][c];
                }
                b[r] = b[r] - factor * b[row];
            }
        }

        // Back substitution
        for (int i = n - 1; i >= 0; --i) {
            // find first non-zero coefficient for pivot column
            int pivcol = -1;
            for (int j = 0; j < n; ++j) {
                if (A[i][j] == fraction(0)) continue; pivcol = j; break;
            }
            // Normally pivcol == i for well-posed systems
            fraction sum = fraction(0);
            for (int j = pivcol + 1; j < n; ++j) {
                sum = sum + A[i][j] * x[j];
            }
            if (A[i][pivcol] == fraction(0)) {
                // free variable; set to 0
                x[i] = fraction(0);
            } else {
                x[pivcol] = (b[i] - sum) / A[i][pivcol];
            }
        }
        // For any remaining variables that didn't get set explicitly
        for (int i = 0; i < n; ++i) {
            // If x[i] is unset, it is default-initialized as 0 by caller
        }
    }

public:
    // 设置电阻网络，构建矩阵A和C。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;

        // Build adjacency (m x n) and conduction (m x m)
        adjacency = matrix(connection_size, interface_size);
        conduction = matrix(connection_size, connection_size);
        for (int e = 0; e < connection_size; ++e) {
            int u = from[e] - 1; // 0-based
            int v = to[e] - 1;
            // adjacency row e: +1 at u, -1 at v
            adjacency(e + 1, u) = fraction(1);
            adjacency(e + 1, v) = fraction(0) - fraction(1);
            // conduction diagonal g = 1 / R
            fraction g = fraction(1) / resistance[e];
            conduction(e + 1, e) = g;
        }
    }

    ~resistive_network() = default;

    // 等效电阻
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1 - 1;
        int b = interface_id2 - 1;
        int n = interface_size;
        if (!(0 <= a && a < n && 0 <= b && b < n && a <= b)) {
            throw resistive_network_error();
        }
        if (n == 1) return fraction(0);
        int k = n - 1;
        // Build reduced conductance matrix G (k x k)
        fraction **G = new fraction *[k];
        for (int i = 0; i < k; ++i) G[i] = new fraction[k];
        build_reduced_conductance(G);
        // Build current injection vector I (k)
        fraction *I = new fraction[k];
        for (int i = 0; i < k; ++i) I[i] = fraction(0);
        if (a != n - 1) I[a] = I[a] + fraction(1);
        if (b != n - 1) I[b] = I[b] - fraction(1);
        // Solve G * U = I
        fraction **A = new fraction *[k];
        for (int i = 0; i < k; ++i) {
            A[i] = new fraction[k];
            for (int j = 0; j < k; ++j) A[i][j] = G[i][j];
        }
        fraction *U = new fraction[k];
        for (int i = 0; i < k; ++i) U[i] = fraction(0);
        solve_linear_system(k, A, I, U);
        fraction Ua = (a == n - 1) ? fraction(0) : U[a];
        fraction Ub = (b == n - 1) ? fraction(0) : U[b];
        fraction req = Ua - Ub;
        // Cleanup
        for (int i = 0; i < k; ++i) { delete [] G[i]; delete [] A[i]; }
        delete [] G; delete [] A; delete [] I; delete [] U;
        return req;
    }

    // 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size的电压为0。
    fraction get_voltage(int id, fraction current[]) {
        int n = interface_size;
        if (!(1 <= id && id < n)) throw resistive_network_error();
        int k = n - 1;
        // Build reduced conductance matrix G (k x k)
        fraction **G = new fraction *[k];
        for (int i = 0; i < k; ++i) G[i] = new fraction[k];
        build_reduced_conductance(G);
        // RHS is current for nodes 1..n-1
        fraction *I = new fraction[k];
        for (int i = 0; i < k; ++i) I[i] = current[i];
        // Solve
        fraction **A = new fraction *[k];
        for (int i = 0; i < k; ++i) {
            A[i] = new fraction[k];
            for (int j = 0; j < k; ++j) A[i][j] = G[i][j];
        }
        fraction *U = new fraction[k];
        for (int i = 0; i < k; ++i) U[i] = fraction(0);
        solve_linear_system(k, A, I, U);
        fraction ans = U[id - 1];
        for (int i = 0; i < k; ++i) { delete [] G[i]; delete [] A[i]; }
        delete [] G; delete [] A; delete [] I; delete [] U;
        return ans;
    }

    // 在给定节点电压U的前提下，返回电阻网络的功率。
    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (int e = 0; e < connection_size; ++e) {
            int u = -1, v = -1;
            for (int j = 0; j < interface_size; ++j) {
                fraction val = adjacency(e + 1, j);
                if (val == fraction(1)) u = j;
                else if (val == fraction(0) - fraction(1)) v = j;
            }
            fraction g = conduction(e + 1, e);
            fraction diff = voltage[u] - voltage[v];
            total = total + g * diff * diff;
        }
        return total;
    }
};


#endif // SRC_HPP

