#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iostream>


using namespace std;

const double PROB_SMOOTH = 1e-7;
const double GLOBALProbabilityForEmpty = 0.4;
const double HMMAlignmentModelSmoothFactor = 0.2;


void
readfile_ch(string str, vector<vector<int>> &vv, unsigned int &cnt, unsigned int &max_line, string output_vcb_file) {

    ifstream fin(str);
    string line;
    string word;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int cnt_line = 0;
    max_line = 0;
    unordered_map<string, int> word2int;

    word2int.insert({"null_wangql", j});//先增加一个null_wangql,主动加空值的一方其实是目标语言（噪声信道）
    ++j;
    ++j;//源语言序号第一个是0，然后从2开始

    while (getline(fin, line)) {
        istringstream linestream(line);
        vv.push_back(vector<int>());

        cnt_line = 0;
        vv[i].push_back(word2int["null_wangql"]);
        ++cnt_line;//每行都增加一个null_wangql

        while (linestream >> word) {
            ++cnt_line;
            auto ret = word2int.insert({word, j});
            if (ret.second) {
                ++j;
            }
            vv[i].push_back(word2int[word]);
        }
        if (cnt_line > max_line) {
            max_line = cnt_line;
        }
        ++i;
    }
    cnt = j;//这里不需要减1是因为虽然序号从2算起，但是开头还有一个0抵消了
    fin.close();

    ofstream fout(output_vcb_file);

    for (auto &wi : word2int) {
        fout << wi.second << " " << wi.first << endl;
    }

    fout.close();
}


void
readfile_en(string str, vector<vector<int>> &vv, unsigned int &cnt, unsigned int &max_line, string output_vcb_file) {
    ifstream fin(str);
    string line;
    string word;
    int i = 0;
    int j = 2;//目标语言序号没有0的存在，直接从2开始
    unsigned int cnt_line = 0;
    max_line = 0;//length of max_line
    unordered_map<string, int> word2int;

    while (getline(fin, line)) {
        istringstream linestream(line);
        vv.push_back(vector<int>());
        cnt_line = 0;
        while (linestream >> word) {
            ++cnt_line;
            auto ret = word2int.insert({word, j});
            if (ret.second) {
                ++j;
            }
            vv[i].push_back(word2int[word]);
        }
        if (cnt_line > max_line) {
            max_line = cnt_line;
        }
        ++i;
    }
    cnt = j - 1;//因为序号从2算起，所以个数要减1
    fin.close();

    ofstream fout(output_vcb_file);

    for (auto &wi : word2int) {
        fout << wi.second << " " << wi.first << endl;
    }

    fout.close();
}

void init_t(vector<vector<int>> &ch, vector<vector<int>> &en, unordered_map<int, unordered_map<int, double>> &t,
            unordered_map<int, unordered_map<int, double>> &t_count, unsigned int cnt) {
    for (unsigned int i = 0; i != ch.size(); ++i) {
        for (auto &c : ch[i]) {
            if (t.find(c) == t.end()) {
                t[c] = unordered_map<int, double>();
                t_count[c] = unordered_map<int, double>();
            }
            for (auto &e : en[i]) {
                t[c][e] = 1.0 / cnt;
                t_count[c][e] = 0.0;
            }
        }
    }
}

void read_t(string str, unordered_map<int, unordered_map<int, double>> &t,
            unordered_map<int, unordered_map<int, double>> &t_count) {

    ifstream fin(str);
    string line;
    string word;
    vector<string> vs(3);
    int ivs = 0;
    while (getline(fin, line)) {
        istringstream linestream(line);
        ivs = 0;
        while (linestream >> word) {
            vs[ivs++] = word;
        }
        if (t.find(stoi(vs[0])) == t.end()) {
            t[stoi(vs[0])] = unordered_map<int, double>();
            t_count[stoi(vs[0])] = unordered_map<int, double>();
        }
        t[stoi(vs[0])][stoi(vs[1])] = stof(vs[2]);
        t_count[stoi(vs[0])][stoi(vs[1])] = 0.0;//t_count结构跟t一样，但是初始化的值为0
    }
    fin.close();

}


void cre_net_n(double **(&net_n), const int &max_line_ch, const int &max_line_en) {//翻译
    int I = 2 * (max_line_ch - 1);
    int J = max_line_en;
    net_n = new double *[I];
    for (int i = 0; i != I; ++i) {
        net_n[i] = new double[J];
    }
    for (int i = 0; i < I; ++i) {
        for (int j = 0; j < J; ++j) {
            net_n[i][j] = 0.0;
        }
    }
}

void del_net_n(double **(&net_n), const int &max_line_ch) {

    int I = 2 * (max_line_ch - 1);

    for (int i = 0; i < I; ++i) {
        delete[] net_n[i];
    }

    delete[] net_n;
}

void ini_net_n(ofstream &fout, double **(&net_n), const int &cnt_line_ch, const int &cnt_line_en,
               unordered_map<int, unordered_map<int, double>> &t, vector<int> &v_ch, vector<int> &v_en) {

    int I = 2 * (cnt_line_ch - 1);
    int J = cnt_line_en;//这里的J是一倍的目标语言长度
    double sum = 0.0;

    for (int j = 0; j < J; ++j) {

        sum = 0;

        for (int i = 0; i < I / 2; ++i) {

            if (t.find(v_ch[i + 1]) == t.end() || t[v_ch[i + 1]].find(v_en[j]) == t[v_ch[i + 1]].end()) {
                net_n[i][j] = PROB_SMOOTH;
            } else {
                net_n[i][j] = max(t[v_ch[i + 1]][v_en[j]], PROB_SMOOTH);
            }

            sum += net_n[i][j];
        }
        for (int i = I / 2; i != I; ++i) {

            if (t.find(v_ch[0]) == t.end() || t[v_ch[0]].find(v_en[j]) == t[v_ch[0]].end()) {
                net_n[i][j] = PROB_SMOOTH;
            } else {
                net_n[i][j] = max(t[v_ch[0]][v_en[j]], PROB_SMOOTH);
            }

            sum += net_n[i][j];
        }
        if (sum) {
            for (int i = 0; i < I; ++i) {
                net_n[i][j] /= sum;
            }
        }

    }
}


void cre_net_e(double **(&net_e), const int &max_line_ch) {//转移
    const int I = 2 * (max_line_ch - 1);
    net_e = new double *[I];
    for (int i = 0; i != I; ++i) {
        net_e[i] = new double[I];
    }
    for (int i = 0; i != I; ++i) {
        for (int j = 0; j != I; ++j) {
            net_e[i][j] = 0.0;
        }
    }
}

void del_net_e(double **(&net_e), const int &max_line_ch) {
    const int I = 2 * (max_line_ch - 1);
    for (int i = 0; i != I; ++i) {
        delete[] net_e[i];
    }
    delete net_e;
}

void
ini_net_e(ofstream &fout, double **(&net_e), const int &cnt_line_ch, unordered_map<int, double> &AlCount, bool doInit) {

    int l = cnt_line_ch - 1;
    int I = 2 * l;

    if (doInit) {//第一次迭代

        for (int i = 0; i != I; ++i) {
            for (int j = 0; j != I / 2; ++j) {
                net_e[i][j] = 1.0 / cnt_line_ch;
            }
            for (int j = I / 2; j != I; ++j) {
                if (i % (cnt_line_ch - 1) == j % (cnt_line_ch - 1)) {
                    net_e[i][j] = 1.0 / cnt_line_ch;
                } else {
                    net_e[i][j] = 0.0;
                }
            }
        }

        //归一化
        for (int i = 0; i != I; ++i) {
            double sum = 0.0;
            for (int j = 0; j != I; ++j) {
                sum += net_e[i][j];
            }
            for (int j = 0; j != I; ++j) {
                net_e[i][j] /= sum;
            }
        }

    } else {

        vector<double> al(l, 0.0);
        double sum = 0.0;
        for (int i = 0; i != I; ++i) {
            sum = 0.0;
            for (int j = 0; j != l; ++j) {
                al[j] = AlCount[i % l - j];
                //                fout<<l<<" "<<i<<" "<<j<<" "<<al[j]<<endl;
                sum += al[j];
            }
            for (int j = 0; j != l; ++j) {
                al[j] /= sum;
                //                fout<<l<<" "<<i<<" "<<j<<" "<<al[j]<<endl;
            }

            double p = HMMAlignmentModelSmoothFactor;
            double pp = p / l;
            for (int j = 0; j != l; ++j) {
                al[j] = (1.0 - p) * al[j] + pp;
                //                fout<<l<<" "<<i<<" "<<j<<" "<<al[j]<<endl;
            }

            double sum = 0.0;
            for (int j = 0; j != I; ++j) {
                net_e[i][j] = al[j % l];

                if (j >= l) {
                    if ((i % l) != (j % l)) {
                        net_e[i][j] = 0;
                    } else {
                        net_e[i][j] = GLOBALProbabilityForEmpty;
                    }
                }
                sum += net_e[i][j];
            }


            if (sum) {
                for (int j = 0; j != I; ++j) {
                    net_e[i][j] /= sum;
                }
            } else {
//                cout<<"这不太可能！"<<endl;
                for (int j = 0; j != I; ++j) {
                    net_e[i][j] = 1.0 / l;
                }
            }

        }
    }
}


void makeHMMNetwork(ofstream &fout, vector<int> &v_ch, vector<int> &v_en, double **(&net_n),
                    unordered_map<int, unordered_map<int, double>> &t, double **(&net_e),
                    unordered_map<int, double> &AlCount, bool doInit) {

    int cnt_line_ch = (int) v_ch.size();
    int cnt_line_en = (int) v_en.size();


    ini_net_e(fout, net_e, cnt_line_ch, AlCount, doInit);
    ini_net_n(fout, net_n, cnt_line_ch, cnt_line_en, t, v_ch, v_en);

    //from now on...

}


void cre_AlCount(unsigned int &max_line, unordered_map<int, double> &AlCount) {
    AlCount[0] = 0.0;
    for (int i = 1; i != max_line; ++i) {
        AlCount[i] = 0.0;
        AlCount[-i] = 0.0;
    }
}


void cre_ai(unsigned int &max_line, unordered_map<int, vector<double>> &ai) {

    for (int i = 1; i != max_line; ++i) {
        ai[i] = vector<double>(2 * i, 0.0);
    }
}

void reset_ai_count(unsigned int &max_line, unordered_map<int, vector<double>> &ai_count) {//置0

    for (int i = 1; i != max_line; ++i) {
        ai_count[i].assign(ai_count[i].size(), 0.0);
    }
}


void addAlCount(ofstream &fout, const vector<int> &v_ch, const vector<int> &v_en, const vector<double> &E,
                unordered_map<int, double> &AlCount) {

    const int l = (int) v_ch.size() - 1;
    const int I = 2 * l;
    auto ep = E.begin();
    double mult = 1.0;
    mult *= l;
    int i_bef_real = 0;

    double np0c = 0.0;

    for (int i = 0; i != l; ++i) {//只算前半部分,但是在这里只是写个l是错误的，因为这样的话ep加的就错了
        for (int i_bef = 0; i_bef != I; ++i_bef, ++ep) {
            i_bef_real = (i_bef % l);
            //            AlCount[l][i_bef_real-i+l-1]+= *ep * mult;//原来的下标是关于y轴对称的，加上l-1后都位于y轴右边，才能用 vector 的坐标访问
            AlCount[i_bef_real - i] += *ep * mult;
            np0c += *ep * mult;
        }
    }

}


void addAiBi(ofstream &fout, vector<int> &v_ch, vector<double> &gamma, unordered_map<int, vector<double>> &ai,
             unordered_map<int, vector<double>> &bi) {
    const int I = 2 * ((int) v_ch.size() - 1);
    auto gp1 = gamma.begin();
    auto gp2 = gamma.end() - I;
    for (int i = 0; i != I; ++i, ++gp1, ++gp2) {
        ai[I / 2][i] += *gp1;
        bi[I / 2][i] += *gp2;
    }
}


void ForwardBackwordTraining(ofstream &fout, vector<int> &v_ch, vector<int> &v_en, double **(&net_n), double **(&net_e),
                             vector<double> &gamma, vector<double> &E,
                             unordered_map<int, unordered_map<int, double>> &t_count,
                             unordered_map<int, vector<double>> &ai, unordered_map<int, vector<double>> &bi,
                             bool doInit, vector<double> &betainit_global) {
    const int I = 2 * ((int) v_ch.size() - 1);
    const int J = (int) v_en.size();
    const int N = I * J;
    vector<double> alphainit(I, 1.0);//这里只是占坑，真正的初始化在下面
    vector<double> betainit(I, 1.0);
    vector<double> alpha(N, 0);
    vector<double> beta(N, 0);


    if (doInit) {
        //这里的初始化情况仅仅对于第一次有效
        double sum_alphainit = 0.0;
        for (int i = 0; i < I; ++i) {
            alphainit[i] = (i < I / 2) ? 1 : (2.0 / I);//第一遍初始化分为前半段和后半段
            sum_alphainit += alphainit[i];
        }
        for (int i = 0; i < I; ++i) {
            alphainit[i] /= sum_alphainit;
        }

        double sum_betainit = 0.0;
        for (int i = 0; i < I; ++i) {
            sum_betainit += betainit[i];
        }
        for (int i = 0; i != I; ++i) {
            betainit[i] /= sum_betainit;
        }
        transform(betainit.begin(), betainit.end(), betainit.begin(), bind1st(multiplies<double>(), I));


    } else {
        double sum_alphainit = 0.0;
        for (int i = 0; i != I / 2 + 1; ++i) {//only first empty word can be chosen
            alphainit[i] = ai[I / 2][i];
            sum_alphainit += alphainit[i];
        }

        for (int i = I / 2 + 1; i != I; ++i) {//后半段除了第一个空其余都为0
            alphainit[i] = 0.0;
        }

        for (int i = 0; i != I; ++i) {
            alphainit[i] /= sum_alphainit;
        }

        double sum_betainit = 0.0;
        for (int i = 0; i != I; ++i) {
            betainit[i] = bi[I / 2][i];
            sum_betainit += betainit[i];
        }
        for (int i = 0; i != I; ++i) {
            betainit[i] /= sum_betainit;
        }
        transform(betainit.begin(), betainit.end(), betainit.begin(), bind1st(multiplies<double>(), I));
        //这里的目的是让每个元素乘以I，这样所有元素的和为I（2*l）

    }


    for (int i = 0; i < I; ++i) {
        beta[N - I + i] = betainit[i];
        betainit_global[i] = betainit[i];
    }

    int NN = N - I - 1;

    for (int j = J - 2; j >= 0; --j) {
        for (int ti = I - 1; ti >= 0; --ti, --NN) {
            auto next_beta = beta.begin() + (j + 1) * I;//后一列beta值最上面一个
            for (int ni = 0; ni < I; ++ni) {
                beta[NN] += (*next_beta++) * (net_e[ti][ni]) * (net_n[ni][j + 1]);
            }
        }
    }


    for (int i = 0; i < I; ++i) {
        alpha[i] = alphainit[i] * net_n[i][0];
    }

    auto cur_alpha = alpha.begin() + I;

    auto cur_beta = beta.begin() + I;

    E.resize(I * I);
    fill(E.begin(), E.end(), 0.0);

    for (int j = 1; j < J; ++j) {

        for (int ti = 0; ti < I; ++ti, ++cur_alpha, ++cur_beta) {
            auto prev_alpha = alpha.begin() + I * (j - 1);//迭代器
            auto this_node = net_n[ti][j];//double类型的翻译概率
            for (int pi = 0; pi < I; ++pi, ++prev_alpha) {
                auto alpha_increment = *prev_alpha * net_e[pi][ti] * this_node;
                (*cur_alpha) += alpha_increment;
                E[I * ti + pi] += alpha_increment * (*cur_beta);
            }
        }
    }

    gamma.resize(N);//这个resize很管用
    transform(alpha.begin(), alpha.end(), beta.begin(), gamma.begin(), multiplies<double>());//这个gamma是针对一句话的

    auto ge = gamma.end();
    for (auto gp = gamma.begin(); gp != ge; gp += I) {
        double sum = 0;
        for (auto gval = gp; gval != gp + I; ++gval) {
            sum += *gval;
        }
        if (sum) {
            for (auto gval = gp; gval != gp + I; ++gval) {
                *gval /= sum;
            }
        } else {
            fill(gp, gp + I, 1.0 / I);
        }
    }

    auto gp = gamma.begin();
    for (int j = 0; j != J; ++j) {
        for (int i = 0; i != I; ++i, ++gp) {
            if (*gp > PROB_SMOOTH) {
                if (i >= I / 2) {//i>=I/2的部分表示空
                    t_count[v_ch[0]][v_en[j]] += *gp;
                } else {
                    t_count[v_ch[i + 1]][v_en[j]] += *gp;//v_ch下标加1用来跳过开头的空
                }
            }
        }
    }

    double Esum = 0.0;

    for (auto Ep = E.begin(); Ep != E.end(); ++Ep) {
        Esum += *Ep;
    }
    if (Esum) {
        for (auto Ep = E.begin(); Ep != E.end(); ++Ep) {
            *Ep /= Esum;
            //            fout << "e：" << I << ":" << J << "  " << kk++ << "  " << *Ep << endl;
        }
    } else {//这一句应该不会执行
        for (auto Ep = E.begin(); Ep != E.end(); ++Ep) {
            *Ep = 1.0 / (I * I);
            //            fout << "e：" << I << ":" << J << "  " << kk++ << "  " << *Ep << endl;
//            cout<<"执行了我不认为会执行的一句！"<<endl;
        }
    }


}

void HMMRealViterbi(unordered_map<int, int> &word_freq_ch, unordered_map<int, int> &word_freq_en, ofstream &fout,
                    vector<int> &viterbi_alignment, vector<int> &vitar, vector<int> &v_ch, vector<int> &v_en,
                    double **(&net_n), double **(&net_e), unordered_map<int, vector<double>> &ai, bool doInit,
                    vector<double> &betainit_global) {
    const int l = (int) v_ch.size() - 1;
    const int I = 2 * l;
    const int J = (int) v_en.size();
    const int N = I * J;
    vector<double> alpha(N, -1);
    vector<double *> bp(N, (double *) 0);//用以存储当前节点取最大值时前一节点
    vitar.resize(J);

    vector<double> alphainit(I, 1.0);//这里只是占坑，真正的初始化在下面


    if (doInit) {
        //这里的初始化情况仅仅对于第一次有效
        double sum_alphainit = 0.0;
        for (int i = 0; i < I; ++i) {
            alphainit[i] = (i < I / 2) ? 1 : (2.0 / I);//第一遍初始化分为前半段和后半段
            sum_alphainit += alphainit[i];
        }
        for (int i = 0; i < I; ++i) {
            alphainit[i] /= sum_alphainit;
        }
    } else {
        double sum_alphainit = 0.0;
        for (int i = 0; i != I; ++i) {
            alphainit[i] = ai[I / 2][i];
            sum_alphainit += alphainit[i];
        }
        for (int i = 0; i != I; ++i) {
            alphainit[i] /= sum_alphainit;
        }
    }


    for (int i = 0; i < I; ++i) {
        alpha[i] = alphainit[i] * net_n[i][0];
        //        cout << I << ":" << J << "   " << i << "    beta    " << beta[i] <<"   betainit    "<<betainit[i]<<"    net_n    "<<net_n[i][0]<< endl;

        if (i > I / 2) {//only first empty word can be chosen
            alpha[i] = 0;
        }
//
//        fout<<alphainit[i]<<"*"<<net_n[i][0]<<endl;//由此发现查找到的不同之处出现在net_n[i][0]上面
//        fout<<I<<"   alpha["<<i<<"]   "<<alpha[i]<<endl;
        bp[i] = 0;
    }


    auto cur_alpha = alpha.begin() + I;
    //    auto cur_bp=bp.begin()+I;
    double **cur_bp = (&*bp.begin()) + I;
    for (int j = 1; j < J; ++j) {
        for (int ti = 0; ti < I; ++ti, ++cur_alpha, ++cur_bp) {
            double *prev_alpha = &*(alpha.begin() + I * (j - 1));

            //            auto prev_alpha=alpha.begin()+I*(j-1);
            double this_node = net_n[ti][j];//翻译概率

//            if(word_freq_ch[v_ch[ti]]>2 || ti >=l){//加了一条限制,但是我这个限制好像并没有起作用

            for (int pi = 0; pi < I; ++pi, ++prev_alpha) {
                //                cout<<*prev_alpha<<endl;
                //                double test=*prev_alpha;
                double alpha_increment =
                        *prev_alpha * net_e[pi][ti] * this_node;//之前在这里访问错误是因为上面的pi最后一部分忘记写++pi了，因此这里不能少
                if (alpha_increment > *cur_alpha) {
                    (*cur_alpha) = alpha_increment;
                    (*cur_bp) = prev_alpha;//存放是指针
                }
            }

//            }
        }
    }





//    vector<double> betainit(I, 1.0);
//
//
//    for (int i=0; i<I; ++i) {
//        alpha[N-I+i]*=betainit[i];
//    }
    //上面这一段没有任何意义，但是对于原始HMM算法,最后一列确实没有乘betainit
    //对于下面乘以betainit的情况，可以认为乘的是那个隐变量出现的概率

    for (int i = 0; i < I; ++i) {
        alpha[N - I + i] *= betainit_global[i];
        //cout<<betainit_global[i]<<endl;
    }


    int j = J - 1;
    cur_alpha = alpha.begin() + j * I;
    vitar[J - 1] = int(max_element(cur_alpha, cur_alpha + I) - cur_alpha);//max_element返回的是迭代器，不是元素，这里返回的是在这一列的索引
    while (bp[vitar[j] + j * I]) {//bp里面放的是此点的前一个点的指针，vitar[j]+j*I得到的是全局的索引
        cur_alpha -= I;
        vitar[j - 1] = int(bp[vitar[j] + j * I] - (&*cur_alpha));//这样减掉本列的初始值确实得到本列的偏移值
        --j;
    }

    viterbi_alignment.resize(J + 1);//注意J+1个长度
    for (int j = 1; j <= J; ++j) {
        viterbi_alignment[j] = vitar[j - 1] + 1;
        if (viterbi_alignment[j] > l) {//对的位置大于实际长度，就是对空了
            viterbi_alignment[j] = 0;
        }
    }


    bool firstword = true;
    for (int j = 1; j <= J; ++j) {
        if (viterbi_alignment[j]) {
            if (firstword) {
                firstword = false;
            } else {
                fout << " ";
            }
            fout << viterbi_alignment[j] - 1 << "-" << j - 1;
        }
    }
    fout << endl;

}

void normalizeTableMy(unordered_map<int, unordered_map<int, double>> &t_count) {
    for (auto it1 = t_count.begin(); it1 != t_count.end(); ++it1) {
        double tmp_sum = 0.0;
        for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            tmp_sum += it2->second;
        }
        if (tmp_sum > 0.0) {
            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
                it2->second /= tmp_sum;
            }
        }
    }
}

void normalizeTable(ofstream &fout, unordered_map<int, unordered_map<int, double>> &t_count, vector<int> &nCh,
                    vector<int> &nEn, const int cnt_ch, const int cnt_en, int it = 2) {
    int tr_cnt_ch = cnt_ch + 1;//加个1号占坑
    int tr_cnt_en = cnt_en + 2;//加个0号和1号占坑

    nCh.resize(tr_cnt_ch);
    for (int i = 0; i != tr_cnt_ch; ++i) {
        nCh[i] = 0;
    }

    nEn.resize(tr_cnt_en);
    for (int i = 0; i != tr_cnt_en; ++i) {
        nEn[i] = 0;
    }

    vector<double> total(tr_cnt_ch, 0.0);//这里+1是因为没有index是1的词，所以坑要多一个
    for (auto iter1 = t_count.begin(); iter1 != t_count.end(); ++iter1) {
        for (auto iter2 = (iter1->second).begin(); iter2 != (iter1->second).end(); ++iter2) {
            total[iter1->first] += iter2->second;
            ++nCh[iter1->first];//每个ch对过几个en
            ++nEn[iter2->first];//每个en对过几个ch
        }
    }
    for (int k = 0; k != tr_cnt_ch; ++k) {
        if (nCh[k]) {
            double probMass = (cnt_en - nCh[k]) * PROB_SMOOTH;//不懂是什么意思
            total[k] += total[k] * probMass / (1 - probMass);
        }
    }

    double p;
    for (auto iter1 = t_count.begin(); iter1 != t_count.end(); ++iter1) {
        for (auto iter2 = (iter1->second).begin(); iter2 != (iter1->second).end(); ++iter2) {
            if (total[iter1->first] > 0.0) {
                p = 1.0 * iter2->second / total[iter1->first];
            } else {
                p = 0.0;
            }
            if (p > PROB_SMOOTH) {
                //                if (it>0) {
                //                    iter2->second=p;
                //                }else{
                //                    //这里应该传给t了
                //                }
                iter2->second = p;
            } else {
                iter2->second = 0.0;
            }
        }
    }
    if (it > 0) {
        normalizeTable(fout, t_count, nCh, nEn, cnt_ch, cnt_en, it - 1);
    }
}

void ai_count2ai(unordered_map<int, vector<double>> &ai,
                 unordered_map<int, vector<double>> &ai_count) {//把ai_count中存储的值传递到ai并且清零ai

//    for (auto i1=ai.begin(),i2=ai_count.begin(); i1!=ai.end()&&i2!=ai_count.end(); ++i1,++i2) {
//        for (auto j1=i1->second.begin(),j2=i2->second.begin(); j1!=i1->second.end()&&j2!=i2->second.end(); ++j1, ++j2) {
//            *j1=*j2;
//            *j2=0.0;
//        }
//    }

    for (auto i = ai.begin(); i != ai.end(); ++i) {
        ai[i->first] = ai_count[i->first];
        ai_count[i->first] = vector<double>(2 * (i->first), 0.0);
    }

}

void AlCount_count2AlCount(unordered_map<int, double> &AlCount, unordered_map<int, double> &AlCount_count) {

//    for (auto i1=AlCount.begin(), i2=AlCount_count.begin(); i1!=AlCount.end()&&i2!=AlCount_count.end(); ++i1, ++i2) {
//        i1->second = i2->second;
//        i2->second = 0.0;
//    }

    for (auto i = AlCount.begin(); i != AlCount.end(); ++i) {
        AlCount[i->first] = AlCount_count[i->first];
        AlCount_count[i->first] = 0.0;
    }
}

void unordered_map_t_count2map_t(unordered_map<int, unordered_map<int, double>> &unordered_map_t,
                                 unordered_map<int, unordered_map<int, double>> &unordered_map_t_count) {
//    for (auto i1=unordered_map_t.begin(), i2=unordered_map_t_count.begin(); i1!=unordered_map_t.end()&& i2!=unordered_map_t_count.end(); ++i1, ++i2) {
//        for (auto j1=i1->second.begin(), j2=i2->second.begin(); j1!=i1->second.end()&&j2!=i2->second.end(); ++j1, ++j2) {
//            j1->second = j2->second;
//            j2->second =0.0;
//        }
//    }

    for (auto i = unordered_map_t.begin(); i != unordered_map_t.end(); ++i) {
        for (auto j = i->second.begin(); j != i->second.end(); ++j) {
            j->second = unordered_map_t_count[i->first][j->first];
            unordered_map_t_count[i->first][j->first] = 0.0;
        }
    }
}

void print_t(string output_file, unordered_map<int, unordered_map<int, double>> &unordered_map_t) {
    ofstream fout(output_file);
    for (auto it1 = unordered_map_t.begin(); it1 != unordered_map_t.end(); ++it1) {
        for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            fout << it1->first << " " << it2->first << " " << it2->second << endl;
        }
    }
    fout.close();
}

void PrintAlCount(string output_file, unordered_map<int, double> &AlCount) {
    ofstream fout(output_file);
    for (auto it = AlCount.begin(); it != AlCount.end(); ++it) {
        fout << it->first << " " << it->second << endl;
    }
    fout.close();
}

void Printai(string output_file, unordered_map<int, vector<double>> &ai) {
    ofstream fout(output_file);
    for (auto it1 = ai.begin(); it1 != ai.end(); ++it1) {
        for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            fout << it1->first << " " << *it2 << endl;
        }
    }
    fout.close();
}

void Printbi(string output_file, unordered_map<int, vector<double>> &bi) {
    ofstream fout(output_file);
    for (auto it1 = bi.begin(); it1 != bi.end(); ++it1) {
        for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            fout << it1->first << " " << *it2 << endl;
        }
    }
    fout.close();
}


int main() {

    time_t start, stop;
    vector<vector<int>> vv_ch;
    unsigned int cnt_ch;
    unsigned int max_line_ch;
    readfile_ch("/root/corpus/LDC.nosemi.final.ch", vv_ch, cnt_ch, max_line_ch, "/root/modelh/LDC.nosemi.final.ch.vcb");

    vector<vector<int>> vv_en;
    unsigned int cnt_en;
    unsigned int max_line_en;
    readfile_en("/root/corpus/LDC.final.en", vv_en, cnt_en, max_line_en, "/root/modelh/LDC.final.en.vcb");


    unordered_map<int, unordered_map<int, double>> unordered_map_t, unordered_map_t_count;//t用来存概率，t_count用来存储中间计算用来更新概率的结果
    init_t(vv_ch, vv_en, unordered_map_t, unordered_map_t_count, cnt_en);

    double **net_n;
    cre_net_n(net_n, max_line_ch, max_line_en);

    double **net_e;
    cre_net_e(net_e, max_line_ch);

    vector<double> gamma;//gamma每个句对都不一样，因此需要在每个句对前面进行构造一下
    vector<double> E;//E每个句对都不一样，因此同样需要在每个句对前面进行构造一下
    vector<int> vit;
    vector<int> viterbi_alignment;

    vector<int> nCh;//每个ch对应en的个数
    vector<int> nEn;

    unordered_map<int, double> AlCount, AlCount_count;//用来保存net_e变化的中间结果
    cre_AlCount(max_line_ch, AlCount);
    cre_AlCount(max_line_ch, AlCount_count);//存放中间结果，类似unordered_map_t_count的作用

    unordered_map<int, vector<double>> ai, bi, ai_count, bi_count;//分别用来存放alpha，beta初始化
    cre_ai(max_line_ch, ai);
    cre_ai(max_line_ch, bi);
    cre_ai(max_line_ch, ai_count);//存放中间结果，类似unordered_map_t_count的作用
    cre_ai(max_line_ch, bi_count);

    ofstream fout("/root/modelh/LDC.modelh.align");

    start = time(NULL);
    for (int i = 0; i < vv_ch.size(); ++i) {
        makeHMMNetwork(fout, vv_ch[i], vv_en[i], net_n, unordered_map_t, net_e, AlCount, 1);

        const int I = 2 * ((int) vv_ch[i].size() - 1);
        vector<double> betainit_global(I, 0.0);

        ForwardBackwordTraining(fout, vv_ch[i], vv_en[i], net_n, net_e, gamma, E, unordered_map_t_count, ai, bi, 1,
                                betainit_global);//因为下面没有用到gamma，E,所以这些结果或许仅仅用于更新第二遍以后的概率，而下面的HMMRealViterbi
        //HMMRealViterbi(fout, viterbi_alignment, vit, vv_ch[i], vv_en[i], net_n, net_e, ai, 1, betainit_global);
        addAlCount(fout, vv_ch[i], vv_en[i], E, AlCount_count);
        addAiBi(fout, vv_ch[i], gamma, ai_count, bi_count);

    }//


    normalizeTable(fout, unordered_map_t_count, nCh, nEn, cnt_ch, cnt_en);
//    normalizeTableMy(unordered_map_t_count);

    stop = time(NULL);
    cout << stop - start << endl;

    unordered_map<int, int> word_freq_ch;
    unordered_map<int, int> word_freq_en;

    for (int loop = 0; loop != 4; ++loop) {
        unordered_map_t_count2map_t(unordered_map_t, unordered_map_t_count);
        AlCount_count2AlCount(AlCount, AlCount_count);
        ai_count2ai(ai, ai_count);
        ai_count2ai(bi, bi_count);


        if (3 == loop) {
            print_t("/root/modelh/LDC.final.modelh.t", unordered_map_t);
            PrintAlCount("/root/modelh/LDC.final.modelh.AlCount", AlCount);
            Printai("/root/modelh/LDC.final.modelh.ai", ai);
            Printbi("/root/modelh/LDC.final.modelh.bi", bi);
        }


        for (int i = 0; i < vv_ch.size(); ++i) {
            makeHMMNetwork(fout, vv_ch[i], vv_en[i], net_n, unordered_map_t, net_e, AlCount, 0);
            const int I = 2 * ((int) vv_ch[i].size() - 1);
            vector<double> betainit_global(I, 0.0);
            ForwardBackwordTraining(fout, vv_ch[i], vv_en[i], net_n, net_e, gamma, E, unordered_map_t_count, ai, bi, 0,
                                    betainit_global);
            if (loop == 3) {


                HMMRealViterbi(word_freq_ch, word_freq_en, fout, viterbi_alignment, vit, vv_ch[i], vv_en[i], net_n,
                               net_e, ai, 0, betainit_global);
            }
            addAlCount(fout, vv_ch[i], vv_en[i], E, AlCount_count);
            addAiBi(fout, vv_ch[i], gamma, ai_count, bi_count);
        }
        normalizeTable(fout, unordered_map_t_count, nCh, nEn, cnt_ch, cnt_en);
//        normalizeTableMy(unordered_map_t_count);
        stop = time(NULL);
        cout << stop - start << endl;

    }
    fout.close();

    stop = time(NULL);
    printf("hmm corpus use time:%ld\n", (stop - start));

    return 0;
}
