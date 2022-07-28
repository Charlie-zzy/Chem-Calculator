#include <iostream>
#include <vector>

typedef long long ll;
typedef long double ld;
#define MAXN 256
using namespace std;

// Utils
ll gcd(ll a, ll b)
{
    while (b)
    {
        ll t = a % b;
        a = b, b = t;
    }
    return a;
}
ll lcm(ll a, ll b) { return (a * b) / gcd(a, b); }

// 高精分数实现
class fraction
{
public:
    ll a, b;
    fraction(ll n = 0) { a = n, b = 1; }
    fraction(ll _a, ll _b)
    {
        a = _a, b = _b;
        ll t = gcd(_a, _b);
        if (t)
            a /= t, b /= t;
    }
    fraction operator+=(fraction n)
    {
        fraction t(a * n.b + b * n.a, b * n.b);
        a = t.a, b = t.b;
    }
    fraction operator-=(fraction n)
    {
        fraction t(a * n.b - b * n.a, b * n.b);
        a = t.a, b = t.b;
    }
    fraction operator*=(fraction n)
    {
        fraction t(a * n.a, b * n.b);
        a = t.a, b = t.b;
    }
    fraction operator/=(fraction n)
    {
        fraction t(a * n.b, b * n.a);
        a = t.a, b = t.b;
    }
    fraction operator+(fraction n) { return fraction(a * n.b + b * n.a, b * n.b); }
    fraction operator-(fraction n) { return fraction(a * n.b - b * n.a, b * n.b); }
    fraction operator*(fraction n) { return fraction(a * n.a, b * n.b); }
    fraction operator/(fraction n) { return fraction(a * n.b, b * n.a); }
    bool operator>=(fraction n) { return a * n.b >= b * n.a; }
    bool operator>(fraction n) { return a * n.b > b * n.a; }
    bool operator<=(fraction n) { return a * n.b <= b * n.a; }
    bool operator<(fraction n) { return a * n.b < b * n.a; }
    bool operator==(fraction n) { return a * n.b == b * n.a; }
    bool operator!=(fraction n) { return a * n.b != b * n.a; }
    operator ld() { return (ld)a / (ld)b; }
} mat[MAXN][MAXN], ans[MAXN], elec[MAXN], e;          // 系数矩阵和增广列
short color[10000];                                   // 字符串上色
ll m, n, N, testIdx = 0, cnt = 1, rev, realAns[MAXN]; // 化学式个数、（调试用）字符串映射、颜色映射、（调试用）测试次数
bool isOrdered[MAXN], debugMode = false;
string strMap[MAXN]; // 元素名映射
vector<string> demo = {"CuSO4 + NaOH = Cu(OH)2 + Na2SO4",
                       "Cu+HNO3 = Cu(NO3)2+NO+H2O",
                       "HaoYe = Hao + Ye",
                       "(C2H4)1000000+O2=CO2+H2O",
                       "ABC=A2+B3+C7",
                       "K2MnO4+H2C2O4+H2SO4 = K2SO4+CO2+H2O+MnO2",
                       "KMnO4+H2C2O4+H2SO4 == K2SO4+CO2+H2O+MnO2",
                       "Ca [2+] + OH [-] = Ca(OH)2",
                       "CH4+O[2+]+OH[-]=CO3[2-]+H2O",
                       "CH4+OH[-]=CO2+H2O",
                       "CO3[2-]+H2O=CH3+O[-]",
                       "4NaHCO3+3Ca(OH)2=CaCO3+Ca(HCO3)2+H2O+NaOH"};

ll toNumber(string s, ll def = 0)
{
    if (s == "")
        return def;
    ll r = (s[0] == '-' ? -1 : s[0] - '0');
    for (ll i = 1; i < s.length(); i++)
        r = r * 10 + ll(s[i] - '0');
    return r;
}
void sayIn(bool newline = false)
{
    if (newline)
        cnt++;
    cout << "\033[33m\033[1mIn [\033[0m\033[36m" << cnt << "\033[1m\033[33m]" << (cnt < 10 ? " " : "") << "\033[0m: ";
}
void sayOut()
{
    cout << "\033[33m\033[1mOut[\033[0m\033[36m" << cnt << "\033[1m\033[33m]" << (cnt < 10 ? " " : "") << "\033[0m: ";
    cnt++;
}
bool isIn(char c, string pattern) { return pattern.find(c) != pattern.npos; }
bool isUpperCase(char c) { return 'A' <= c && c <= 'Z'; }
bool isLowerCase(char c) { return 'a' <= c && c <= 'z'; }
bool isNumber(char c) { return '0' <= c && c <= '9'; }
bool isValid(char c) { return isUpperCase(c) || isLowerCase(c) || isNumber(c) || isIn(c, "[-]=(+)"); }

enum // 有逼格的颜色代码
{
    PLUS = -1,
    EQUAL = -2,
    BRACE = -3,
    ANS = -4,
    NUM = -5,
    DELETE = -9,
    XX = -10
};
enum // 有逼格的返回值
{
    OK,
    ERROR,
    ERROR_INFINITY_SOLUTION,
    ERROR_NO_SOLUTION,
    ERROR_DIVIDE_BY_ZERO
};

string clr(ll c) // 颜色映射
{
    switch (c)
    {
    case ANS:
        return "\033[1m\033[33m";
    case DELETE:
        return "\033[9m";
    default:
        return "\033[0m";
    }
}
string formatNum(char c, bool up)
{
    switch (c)
    {
    case '0':
        return (up ? "⁰" : "₀");
    case '1':
        return (up ? "¹" : "₁");
    case '2':
        return (up ? "²" : "₂");
    case '3':
        return (up ? "³" : "₃");
    case '4':
        return (up ? "⁴" : "₄");
    case '5':
        return (up ? "⁵" : "₅");
    case '6':
        return (up ? "⁶" : "₆");
    case '7':
        return (up ? "⁷" : "₇");
    case '8':
        return (up ? "⁸" : "₈");
    case '9':
        return (up ? "⁹" : "₉");
    case '+':
        return (up ? "⁺" : "₊");
    case '-':
        return (up ? "⁻" : "₋");
    }
}
ll getStringIndex(string s) // 元素映射
{
    for (ll i = 0; i < m; i++)
        if (s == strMap[i])
            return i;
    strMap[m++] = s;
    return m - 1;
}
void clear() // 暴力清零
{
    for (ll i = 0; i < MAXN; i++)
    {
        for (ll j = 0; j < MAXN; j++)
            mat[i][j] = 0;
        ans[i] = realAns[i] = color[i] = elec[i] = 0;
        isOrdered[i] = false;
        strMap[i] = "";
    }
    m = n = N = e = 0;
}

// Main Function
string input(ll &cnt)
{
    while (true)
    {
        sayIn();
        string s = "", t; // cin -> t ~> s -> return
        int status = OK;  // 合法判断
        bool eq = false;  // 是否已出现 '='
        if (testIdx)
        {
            t = demo[demo.size() - testIdx];
            cout << t << endl;
            testIdx--;
        }
        else
            getline(cin, t);

        if (t == "help" || t == "Help" || t == "HELP" || t == "?" || t == "/?")
        {
            cout << endl
                 << clr(ANS) << "Help" << clr(XX)
                 << ": You can type in some "
                 << "\033[33mChemical Formulas" << clr(XX)
                 << " and press Enter." << endl
                 << "      Then I'll automatically balance it for you.\n\n"
                 << clr(ANS) << "e.g." << clr(XX) << ": ";
            for (ll i = 0; i < demo.size(); i++)
                cout << demo[i] << "\n      ";
            cout << endl;
            cnt++;
            continue;
        }

        if (t == "debug" || t == "DEBUG")
        {
            debugMode ^= 1;
            cout << (debugMode ? "\n\033[32mDEBUG ON\033[0m\n\n" : "\n\033[31mDEBUG OFF\033[0m\n\n");
            cnt++;
            continue;
        }

        if (t == "test")
        {
            testIdx = demo.size();
            cout << endl;
            cnt++;
            continue;
        }

        for (ll i = 0; i < t.length() && status == OK; i++)
        {
            if (t[i] == ' ') // 忽略空格
                continue;
            if (t[i] == '=')
            {
                if (eq) // 出现多个 '='
                {
                    if (s.back() != '=') // awa = awa = awa [不合法]
                        status = ERROR;
                    continue; // awa == awa [合法]
                }
                else // 刚出现 '='
                    eq = true;
            }
            if (s == "+" || s == "=" || !isValid(t[i])) // 字符合法判定
                status = ERROR;
            s.push_back(t[i]);
        }
        if (status == OK && eq)
            return s + '+'; // 防止漏尾部字符

        // on ERROR
        if (t != "")
            cout << "Invalid Input. Type 'help' or '?' to get help.\n\n";
    }
}
int subProcess(string s, const ll k0, const ll offset)
{
    if (!k0) // 记录反转符号位置
        rev = N;
    ll len = s.length(), L, R, ptr, k1, k2 = 1;

    // 获取系数
    for (ptr = 0; isNumber(s[ptr]); ptr++)
        ;
    k1 = toNumber(s.substr(0, ptr), -1);
    if (k1 != -1)
        isOrdered[N] = true, realAns[N] = k1;
    L = ptr;

    // 计算电子
    if (s[len - 1] == ']')
    {
        for (ptr = len - 2; s[ptr] != '['; ptr--)
            color[ptr + offset] = NUM;
        color[len - 1 + offset] = color[ptr + offset] = BRACE;
        if (!isIn(s[len - 2], "+-"))
        {
            cout << "Invalid Input! Use [114+] or [514-] instead of [1919].\n";
            cin.get();
            system("exit");
        }
        elec[N] = k0 * toNumber(s.substr(ptr + 1, len - ptr - 3), 1) * (s[len - 2] == '-' ? -1 : 1);
        R = ptr;
    }
    else
        R = len, elec[N] = 0;

    // 处理式子
    for (ll i = L; i < R; i++)
    {
        if (isUpperCase(s[i]))
        {
            for (ptr = i + 1; isLowerCase(s[ptr]); ptr++)
                ;
            ll idx = getStringIndex(s.substr(i, ptr - i));
            for (i = ptr; isNumber(s[ptr]); ptr++)
                color[i + offset] = NUM;
            if (k1 != -1)
                ans[idx] -= k0 * k1 * k2 * toNumber(s.substr(i, ptr - i), 1);
            else
                mat[idx][n] += k0 * k2 * toNumber(s.substr(i, ptr - i), 1);
            i = ptr - 1;
        }
        if (s[i] == '(')
        {
            color[i + offset] = BRACE;
            for (ptr = i + 1; s[ptr] != ')'; ptr++)
                ;
            for (L = ptr + 1; isNumber(s[L]); L++)
                ;
            k2 = toNumber(s.substr(ptr + 1, L - ptr - 1), 1);
        }
        if (s[i] == ')')
        {
            color[i + offset] = BRACE;
            for (k2 = 1, i++; isNumber(s[i]); i++)
                color[i + offset] = NUM;
            i--;
        }
    }
}
int process(string S)
{
    ll len = S.length();
    clear();

    // 预处理：末尾必须有 '+' ，否则不会计算末尾化学式
    if (S[len - 1] != '+')
        S.push_back('+');

    // 预处理：如果全小写则自动转大写
    bool haveCapital = false;
    for (ll i = 0; i < len; i++)
        haveCapital |= isUpperCase(S[i]);
    if (!haveCapital)
        for (ll i = 0; i < len; i++)
            if (isLowerCase(S[i]))
                S[i] += 'A' - 'a';
    ll L = 0, R = 0, k0 = 1;
    for (bool inE = false; R < len; R++)
    {
        if (isIn(S[R], "+=") && !inE)
        {
            color[R] == (S[R] == '+' ? PLUS : EQUAL);
            subProcess(S.substr(L, R - L), k0, L);
            if (!isNumber(S[L]))
                n++;
            N++;
            L = R + 1;
        }
        k0 = (S[R] == '=' ? -1 : k0);
        inE = (S[R] != '[' ? S[R] != ']' ? inE : false : true);
    }
}
int solve() // 高斯列主元消元
{
    ll j;
    for (ll i = 0; i < m; i++)
    {
        for (j = i; j < m; j++) // 找列主元
            if (mat[j][i] != fraction())
            {
                swap(mat[i], mat[j]);
                swap(ans[i], ans[j]);
                break;
            }
        if (j == m) // 全是零
            continue;

        for (j = i + 1; j < m; j++)
        {
            fraction t = mat[j][i] / mat[i][i]; // ji / ii
            for (ll k = i; k < n; k++)
                mat[j][k] -= mat[i][k] * t; // jk -= jk * (ji / ii)
            ans[j] -= ans[i] * t;
        }
    }

    ll uke = 0;
    for (ll i = 0; i < N; i++)
        uke += isOrdered[i];
    if (m + uke < n - 1) // 秩都不够还想让我解方程？
        return ERROR_INFINITY_SOLUTION;

    for (ll i = n; i < m; i++) // 这些行应该全是零
        if (ans[i] != fraction())
            return ERROR_NO_SOLUTION;

    if (mat[n - 1][n - 1] == fraction()) // 化学方程式肯定会多一个不定元的说
        mat[n - 1][n - 1] = ans[n - 1] = 1, m++;

    for (ll i = n - 1; i >= 0; i--)
    {
        for (ll j = i + 1; j < n; j++)
            ans[i] -= ans[j] * mat[i][j];
        if (mat[i][i] == fraction())
            return (ans[i] == fraction() ? ERROR_INFINITY_SOLUTION : ERROR_NO_SOLUTION);
        ans[i] /= mat[i][i];
    }
    return OK;
}
void debugPrint() // （debug）输出矩阵
{
    if (!debugMode)
        return;
    cout << "[DEBUG INFORMATION]\n\n\033[32m charac \033[0m|\t\033[32m";
    for (ll i = 0; i < n; i++)
        cout << i << "\t";
    cout << "\033[0m|\t\033[34mans\033[0m\n";
    for (ll i = 0; i < 8 * n + 34; i++)
        cout << "-";
    cout << endl;
    for (ll j = 0; j < m; j++)
    {
        cout << "\033[32m " << (strMap[j] == "" ? "*" : strMap[j]) << "\033[0m\t|\t";
        for (ll i = 0; i < n; i++)
        {
            if (mat[j][i] == fraction())
                cout << "\033[8m";
            cout << mat[j][i] << "\033[0m\t";
        }
        cout << "|\t\033[34m" << (ans[j] < fraction() ? "-" : "") << abs(ans[j].a);
        if (ans[j].b != 1)
            cout << "/" << abs(ans[j].b);
        cout << "\033[0m\n";
    }
    for (ll i = 0; i < 8 * n + 34; i++)
        cout << "-";
    cout << "\n\n\n   \033[32mN\033[0m    |\033[32m\t";
    for (ll i = 0; i < N; i++)
        cout << i << "\t";

    cout << "\033[0m|\n";
    for (ll i = 0; i < 8 * N + 34; i++)
        cout << "-";
    cout << "\n\033[36m realAns\033[0m|\t";
    for (ll i = 0; i < N; i++)
        cout << realAns[i] << "\033[0m\t";
    cout << "\033[0m|\t\033[36m0\n\033[31m elect  \033[0m|\t\033[31m";
    for (ll i = 0; i < N; i++)
        cout << elec[i] << "\t";
    cout << "\033[0m|\t\033[31m" << e << "\033[0m\n";
    cout << endl;
}
void formatPrint(const string s)
{
    sayOut();
    for (ll i = 0, L = 0, ptr; i < N; i++)
    {
        cout << clr(realAns[i] ? ANS : DELETE) << realAns[i] << " " << clr(XX);
        while (isNumber(s[L]))
            L++;
        ptr = L;
        for (bool inE = false;; ptr++, inE = (s[ptr] != '[' ? s[ptr] != ']' ? inE : false : true))
            if (isIn(s[ptr], "+=") && !inE)
            {
                for (ll j = L; j < ptr; j++)
                {
                    if (s[j] == '[')
                    {
                        inE = true;
                        continue;
                    }
                    if (s[j] == ']')
                        break;
                    cout << clr(realAns[i] ? color[j] : DELETE);
                    if (color[j] == NUM)
                        cout << formatNum(s[j], inE);
                    else
                        cout << s[j];
                    cout << clr(XX);
                }
                if (s[ptr] == '=' && e != fraction())
                {
                    if (e.a < 0)
                        cout << " " << clr(PLUS) << "+ " << clr(ANS) << -e.a << "\033[31me- " << clr(XX) << clr(s[ptr] == '=' ? EQUAL : PLUS) << s[ptr] << clr(XX) << " ";
                    else
                        cout << " " << clr(s[ptr] == '=' ? EQUAL : PLUS) << s[ptr] << clr(XX) << " " << clr(ANS) << e.a << "\033[31me- " << clr(XX) << clr(PLUS) << "+ ";
                }
                else if (ptr != s.length() - 1)
                    cout << " " << clr(s[ptr] == '=' ? EQUAL : PLUS) << s[ptr] << clr(XX) << " ";
                break;
            }
        L = ptr + 1;
    }
}

int main()
{
    ios::sync_with_stdio(false);
    system("title Chemical Calculator v2.4");
    system("chcp 65001");
    system("cls");
    cout << " - Chemical Calculator Made By Charlie - \n"
         << "Command: help, debug, test\n\n";

    while (true)
    {
        string s = input(cnt);
        process(s);

        ll status = solve();
        if (status != OK) // 尝试加上电子再算一遍
        {
            ll elecId = getStringIndex("~elec~");
            for (ll i = 0; i < N; i++)
                if (isOrdered[i])
                    ans[elecId] -= elec[i];
                else
                    mat[elecId][i] = elec[i];
            status = solve();
        }

        ll denominator = 1;
        bool isAllZero = true, divZero = false;
        for (ll i = 0; i < n; i++)
            denominator = lcm(denominator, ans[i].b), isAllZero &= ans[i] == fraction(), divZero |= ans[i].b == 0;
        for (ll i = 0, j = 0; i < n; i++, j++)
        {
            while (isOrdered[j])
                j++;
            realAns[j] = ans[i].a * denominator / ans[i].b;
            if (realAns[i] < 0)
                status = ERROR_INFINITY_SOLUTION;
        }
        for (ll i = 0; i < N; i++)
            e += (i >= rev ? -1 : 1) * elec[i] * realAns[i];
        if (divZero)
            status = ERROR_DIVIDE_BY_ZERO;
        if (isAllZero)
            status = ERROR_NO_SOLUTION;

        if (status == OK || debugMode)
        {
            formatPrint(s);
            debugPrint();
            cout << "\n\n";
        }
        else
            switch (status)
            {
            case ERROR_INFINITY_SOLUTION:
                cout << "\033[31mError\033[0m\n\nThere are infinity solution for this formula.\nMaybe you can try to specific some parameters.\n\ne.g.:  NaHCO3+ Ca(OH)2=CaCO3+Ca(HCO3)2+H2O+NaOH -> \033[31mError\033[0m\n      \033[33m\033[1m4\033[0mNaHCO3+\033[33m\033[1m3\033[0mCa(OH)2=CaCO3+Ca(HCO3)2+H2O+NaOH -> \033[32mOK\033[0m\n\n";
                break;
            case ERROR_NO_SOLUTION:
                cout << "\033[31mError\033[0m\n\nThere are no solution for this formula.\nCheck your input carefully! :P\n\n";
                break;
            case ERROR_DIVIDE_BY_ZERO:
                cout << "\033[31mError\033[0m\n\nA number is divided by zero.\nWhat the fuck??\n\n";
                break;
            }
    }
    return 0;
}
