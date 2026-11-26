# 数式と単位

記号Tは損失エネルギー、Eは電子の初期運動エネルギー、qは運動量移行をhbarで割った波数である。TとEはeV、qはÅ^-1。

## E1: 振動子の損失関数

\[
L(q,T)=\operatorname{Im}\left[-\epsilon(q,T)^{-1}\right]
=\Theta(T-E_g)\sum_i
\frac{A_i\Gamma_i T}
{\left[(E_i+\alpha_i Hq^2)^2-T^2\right]^2+\Gamma_i^2T^2}.
\]

\[
H=\frac{\hbar^2}{2m_e}\simeq3.80998212\ \mathrm{eV\,Å^2}.
\]

出典: [R5]式(3a),(3b)。同等の式が[R3]式(3)と[R4]にもある。

| 入力 | 意味 | 単位 |
| --- | --- | --- |
| energy_eV | E_i、q=0での振動子中心 | eV |
| strength_eV2 | A_i、振動子強度 | eV² |
| width_eV | Gamma_i、減衰パラメータ | eV |
| dispersion | alpha_i、分散係数 | 無次元 |
| gap_eV | E_g、損失のしきい値 | eV |

Lは無次元。`loss_function` は \(T\le E_g\) で \(L=0\) を返す。複素応答は E4 で定義する。光学和則と整合する振動子強度は物性入力に依存する。

## E2: バルクDIIMFP

\[
K_{\rm bulk}(E,T)=\frac{1}{\pi a_0E}
\int_{q_-}^{q_+}\frac{L(q,T)}{q}\,dq,
\qquad
q_\pm=\frac{\sqrt E\pm\sqrt{E-T}}{\sqrt H}.
\]

\[
a_0=0.529177210903\ \mathrm{Å}.
\]

出典: [R4]式(12)、[R2]式(1)、[R3]式(2)。この初版の定義域は0<T<E。均質・等方・非相対論的な誘電応答モデルであり、表面や正孔の補正は含まない。

出力Kの単位はÅ^-1 eV^-1。例えばK dT dxは経路長dx、損失幅dTに対する一次の散乱確率である。有限の全経路に対する損失分布そのものではない。

qの上下限はエネルギー・運動量保存から来る。小さいTでの桁落ちを避けるため、下限は以下の同値式で計算する。

\[
q_- = \frac{T}{\sqrt H(\sqrt E+\sqrt{E-T})}.
\]

測度は dq/q = d(log q) である。実装は振動子ごとに計算し、無分散の項には E3 の解析式を使う。分散がある項では s = αHq² とおく。共鳴から遠い区間では log s、共鳴付近では asinh[(s−(T−Eᵢ))/(Γᵢ/2)] を積分変数にする。数値誤差の推定値は物性モデルの誤差を含まない。運動量区間の対数の幅は、端点の対数の差ではなく \(\log(q_+/q_-)=2\operatorname{asinh}\sqrt{(E_0-T)/T}\) と `log1p` で求め、損失が \(E_0\) に近いときの桁落ちを避ける（段階4f）。

## E3: 独立な解析基準

alpha_i=0ならLはqに依存せず、

\[
K_{\rm bulk}(E,T)=\frac{L(0,T)}{\pi a_0E}
\log\frac{q_+}{q_-}.
\]

E3 は無分散の値と規格化を確かめる基準になる。無分散の項は実装も E3 を使うので、この照合は前因子と対数の実装を確かめる。分散のある積分は、段階4aの多倍長参照、許容誤差を変えた収束試験、I2 の13a項との照合で確かめる。さらに、ギャップなし・単一振動子では、

\[
\int_0^\infty T L(0,T)\,dT=\frac{\pi A}{2}.
\]

これは振動子モデルの一次モーメントの解析値である。実材料の電子密度に対応するf和則を入力が満たすかどうかの検証とは区別する。

## E4. 振動子 ELF からの複素誘電応答

時間依存性を \(e^{-i\omega t}\) とする。\(F=-1/\epsilon\) とおくと、正の損失エネルギーで \(\operatorname{Im}F=L\) となる。固定した \(q\) で \(\Omega_i=E_i+\alpha_i Hq^2\)、しきい値をかける前の振動子和を \(L_0\) とする。\(F(\infty)=-1\)、\(\operatorname{Im}F(q,-T)=-\operatorname{Im}F(q,T)\) を用いれば、逆誘電関数のクラマース・クローニッヒ関係 [R14 式(4),(7)] は次式になる。

\[
\operatorname{Re}F(q,T)=-1+\frac{2}{\pi}\operatorname{PV}
\int_{E_g}^{\infty}\frac{uL_0(q,u)}{u^2-T^2}\,du.
\tag{E4a}
\]

\(u\)、\(T\)、\(E_g\)、\(\Omega_i\)、\(\Gamma_i\) は eV、\(F\) と \(\epsilon\) は無次元。積分では、E1 の振動子で高エネルギー側を無限遠まで外挿する。CLI の出力格子は評価点を指定する。有限帯域の測定点から再構成する場合は、帯域外の ELF とその不確かさを別途定義する必要がある。振動子強度の再規格化は行わない。[R4 §2.2 式(9),(10) と直後の光学和則の記述]

ギャップなしでは、下半平面に極を持ち E1 の虚部を再現する次の解析式を使う。

\[
F_0(q,T)=-1+\sum_i\frac{A_i}{\Omega_i^2-T^2-i\Gamma_iT}.
\tag{E4b}
\]

ギャップがあるときは、振動子ごとに E4a を計算する。中心が端より低い振動子は、\(x=E_g/u\) で残存帯域を直接積分する。\(\Gamma_i\ge2\Omega_i\) かつ \(\Omega_i^2/\Gamma_i\le E_g\) の過減衰振動子にもこの経路を使う。それ以外では E4b の実部から除去帯域の寄与を引く。主値の極や近接する端点では、分子から定数を引いて残りを積分し、定数項の積分を対数で評価する。例えば除去帯域で \(0<T<E_g\)、\(B(u)=uL_0(q,u)\) とおくと、

\[
\operatorname{PV}\int_0^{E_g}\frac{B(u)}{u^2-T^2}\,du
=\int_0^{E_g}\frac{B(u)-B(T)}{u^2-T^2}\,du
+\frac{B(T)}{2T}\log\left|\frac{E_g-T}{E_g+T}\right|.
\tag{E4c}
\]

非ゼロ強度で \(E_g>0\) の場合、\(\operatorname{Re}F\) は \(T\to E_g\) の両側から \(L_0(q,E_g)/\pi\) を係数として対数発散する。複素応答 API は端と、その数値除外幅内にある正の損失を拒否する（[D7](DIELECTRIC_DESIGN.md#ギャップ端と極)）。全強度がゼロなら常に \(F=-1\)、\(\epsilon=1\) となる。

\[
\epsilon(q,T)=-\frac{1}{F(q,T)}.
\tag{E4d}
\]

E4a を上半平面に延長すると、\(\operatorname{Im}F(q,z)\) は \(\operatorname{Re}z\) と同符号になる。虚軸上では \(F(q,iy)\) は \(F(q,0)\) から \(-1\) に単調減少する。したがって \(F(q,0)>0\) のとき、\(\epsilon\) は上半平面に極を持つ。点評価 `dielectric_response` は \(F(q,0)>\tau_0\) を拒否する。絶対値が \(\tau_0\) 以下なら `static_sign` は `marginal` で、T=0 は極として拒否する。正の T では評価を続ける。\(F(q,0)<0\) かつ鋭いギャップを持つ非ゼロ応答では、ギャップ下に \(F=0\) を満たす点が一つある。これは \(\epsilon\) の実軸上の極であり、点評価で F が厳密にゼロなら `ValueError` を返す。実軸の極に付随する分布としての虚部は点ごとの出力には含まれない。極の近傍で逆数を推定誤差の範囲で解像できなければ、数値未収束を返す。[R14 §I]\(\epsilon=-1/F\) の実部と虚部は、成分ごとに精度を判定する。複素数全体の大きさだけでは判定しない（段階4f）。

| 式 | 実装 | 検証 |
| --- | --- | --- |
| E4a, E4c | `negative_inverse_dielectric_function` | しきい値の上下で、残存帯域と無限遠の尾を別に積分した主値計算と比較 |
| E4b | `negative_inverse_dielectric_function` | ギャップなしの単一 Lorentz 振動子、有限 \(q\)、高エネルギー極限 |
| E4d | `dielectric_response` | 解析式との一致、真空極限、実軸と上半平面の極、ELF の符号 |

## E5. 固定深さ・法線出射での XPS 有効断面積

`simonsen1997_fixed_depth_v0` は、ギャップのない材料と、[M2-2 の設計](M2_2_DESIGN.md) §4 の受理条件を満たす硬いギャップのある材料で、表面法線からの出射角 \(\theta=0\)、生成深さ \(a>0\)、損失 \(0<T=\hbar\omega<E_0\) を扱う。非相対論の \(\hbar v=2\sqrt{HE_0}\) を使うと、\(\Omega_0=\omega_\perp=T/(\hbar v)\) となる。正の損失では式(14)–(16)の \(\delta(\omega_\perp)\) の寄与はゼロで、正孔因子には \(i/\omega_\perp\) が残る。[Simonsen ら (1997) 式(13)–(17)](https://doi.org/10.1103/PhysRevB.56.1612)

原著は、t = 0 に電子と正孔が同じ点に突然生じ、正孔は静止して寿命は無限、電子は直線の軌道を一定の速度で進むとする。速度を一定とみなせるのは、損失が E0 に比べて無視できるときである（p. 1613左、式(3)）。表面は鏡面反射の模型で扱う（p. 1613右）。本実装はこれに従い、速度を初期運動エネルギーから非相対論の式 \(\hbar v=2\sqrt{HE_0}\) で一度だけ決める。式(17)の運動量の上下限は \(E_0-T\) を含むが、速度は変えない。計算できる範囲 \(0<T<E_0\) と、数値で照合した範囲は別である。CLI の `run.json` の `applicability` は、格子の \(\max T/E_0\)、これらの仮定、照合した範囲（[VALIDATION.md](VALIDATION.md) の「段階4b」）を記録し、範囲の外では警告を出す。

`xps_fixed_depth.effective_cross_section_fixed_depth` は、式(17)の上下限を `_parallel_limits` で求め、`_term_integrands` で 13a、13b、14、15-T1/T2/T3、16-e/h を \(\ln k_\parallel\) 上で積分する。E4b の評価点は `epsilon_momentum` で明示する。`k_parallel` は \(q=k_\parallel\)、`k_total` は \(q=\sqrt{k_\parallel^2+\Omega_0^2}\) である。分散がある場合、13a が法線で E2 に一致するのは後者である。[式の対応表 §4.2、Q1](M1_EQUATION_MAP.md)

硬いギャップのある材料（M2-2）では、係数を \(F=-1/\epsilon\) で作る（分母は \(1-F\) だけで、ε の極 F = 0 は正則）。節点の実部は残存帯域の部分分数 D3 をベクトル化して求め、虚部は E1 に \(\theta(T-E_g)\) を掛ける。出力は遅延極限の点値で、広がりは加えない。ギャップ下では \(\delta(F-1)\) の寄与がないことを、外向きに丸めた解析的な上界で証明できた損失だけを計算する。そこでは S、R が実数で、正孔電荷1の `total` の被積分関数が点ごとに 0 になる。端の除外、単調性の条件、誤差の伝播と契約は [M2-2 の設計](M2_2_DESIGN.md)、検証の条件は [事前登録](M2_2_PREREG.md) にある。

媒質・正孔・真空の各項はすべて固体内経路長 \(x=a\) 当たりの eV⁻¹ Å⁻¹ で返す。13a の積分 \(k_\parallel\,dk_\parallel/(\Omega_0^2+k_\parallel^2)\) は無次元である。残りの項でも面積次元の核と \(k_\parallel\,dk_\parallel\) の積は無次元で、追加の前因子 \(1/(2x\omega_\perp)\) も無次元となる。単位は共通の \(1/(\pi a_0 E_0)\) で決まる。単位電荷の正孔差分 `hole` は 14、15-T2、16-h の符号付き和で、`total = no_hole + hole_charge * hole` とする。検査条件は [M2 内部検査 I1–I8](M2_FIRST_CASE.md#41-内部の検査参照ソフトなしで自動テストにする) にある。

数値の精度は、8つの項、`no_hole`、単位電荷の `hole`、`total` のそれぞれについて \(\max(\varepsilon_{\rm abs},\varepsilon_{\rm rel}|K|)\) で判定する。表面反射の係数 \(R=S/(2-S)\) と \(C=2S-R\) は、表面共鳴（\(S\approx2\)）の近くで実部が桁落ちする。そこで、実部と虚部の丸め誤差を別々に見積もり、積分点の座標の丸めも含めて各項に伝える。見積もりが要求精度を超えうる点では、与えた倍精度の値から、係数を有理数で、座標を80桁で計算し直す。それでも要求精度を満たせなければ `FixedDepthConvergenceError` を返す（段階4f）。`quadrature_error` は、全項に共通の誤差の上界である。段階4fからは係数の丸めも含むので、それ以前より大きく出る（ELF データ集の材料で中央値は約20倍）。

## E6. 法線出射の深さ平均（外から与える一定の長さ）

`simonsen1997_depth_average_external_v0`（開発版、M3a）は、E5 の固定深さの11出力を原著の式(18)の重みで深さ方向に平均する。[Simonsen ら (1997) 式(18)](https://doi.org/10.1103/PhysRevB.56.1612)

\[
K_{{\rm sc},j}(T;\lambda)=\int_0^\infty w(x)\,K_j(T;a=x)\,dx,\qquad
w(x)=\frac{x}{\lambda^2}e^{-x/\lambda},\qquad \int_0^\infty w\,dx=1
\]

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| \(\lambda\)（`weight_length_A`） | 重みの長さ。外から与える有限の正の値で、全経路で一定。既定値はない | Å |
| \(x=a\) | 固体内の経路長（法線出射なので生成深さに等しい） | Å |
| \(K_j\) | E5 の固定深さの出力（8項、`no_hole`、`hole`、`total`） | eV⁻¹ Å⁻¹ |
| \(K_{{\rm sc},j}\) | 深さ平均した出力（固体内の経路長当たり、真空中の損失を含む） | eV⁻¹ Å⁻¹ |

- **重み:** 一定率のポアソン過程で、経路 \(x\) の間にちょうど1回非弾性散乱する確率 \((x/\lambda)e^{-x/\lambda}\) を、経路長で規格化したもの。重みの積分 1 は経路長の分布の規格化で、損失エネルギー方向の規格化ではない。全出力に同じ重みを使うので、`total = no_hole + hole_charge * hole` が保たれる。
- **仮定:** 均質・半無限の試料、表面法線からの出射角 0°、弾性散乱なし、一様な深さ分布、単一散乱。軌道ごとに変わる長さ（文献の自己無撞着な長さ）、損失窓の面積、式(19)のスペクトルは扱わない。原著の留保（式(18)は REELS の重みの流用で XPS での妥当性は検討を要すること、\(K_{\rm sc}\) は最初の1回の非弾性散乱だけを表し intrinsic の部分を含むこと）は `warnings` に書く。[原文 Simonsen1997 p.1615]
- **深さの積分:** 数値の深さ求積は使わない。\(\theta=0\) では、13a 以外の項の \(aK_j\) の被積分密度が、深さについて \(\{1,\ e^{-sa},\ e^{-\bar sa},\ e^{-2ka}\}\)（\(s=k+ib\)、\(k=k_\parallel\)、\(b=\Omega_0\)）の線形結合になり、係数は深さによらない。そこで式(18)の深さの積分をラプラス変換で閉じ、\(k_\parallel\) の被積分関数に \(p=1/\lambda\) の有理式の重み（W13、N、W14、W15b、W15c、W16e、W16h）を掛けて、E5 と同じ \(\ln k_\parallel\) 上で積分する。13a は深さによらないので固定深さの値と同じになる。\(a\to0\) で \(aK_j\)（13a を除く）は有限の極限（媒質内の項では 0、真空中の 16-e と 16-h では 0 でない値）を持つので、\(w(x)K_j\) は \(x=0\) で有限で、深さの積分に特異性はない。重みとその誤差の上界は [M3a の設計](M3A_DESIGN.md) §2.2、§3.1 にある。
- **ギャップのある材料:** E5 と同じ受理条件と事前検査（深さによらない）を使う。ギャップ下で正孔電荷が 1 の `total` は、8項の値の和ではなく、節点ごとに8項を足した直接和を積分して作り、その数値誤差を別に評価する（理論値は 0）。
- **精度:** 11出力のそれぞれで、推定誤差が \(\max(\varepsilon_{\rm abs},\varepsilon_{\rm rel}|K_{{\rm sc},j}|)\) 以下になるまでパネルを細分する。満たせなければ `DepthAverageConvergenceError`（`reason` は `budget_not_met`、`panel_floor_exceeds_budget`、`direct_sum_floor_exceeds_budget`、`weight_out_of_range` ほか）で止まる。`quadrature_error` は11出力の推定の最大である。既定は API で \(\varepsilon_{\rm abs}=10^{-12}\) Å⁻¹ eV⁻¹、\(\varepsilon_{\rm rel}=10^{-9}\)、CLI では設定ファイルの `integration`（省略時は \(10^{-9}\)、\(10^{-7}\)）で、limit はともに 200。

実装は `xps_depth_averaged.effective_cross_section_depth_averaged`（深さの閉形式は `_depth_kernel`）、CLI は `xps-depth-averaged`。検証の条件と許容幅は [M3a の事前登録](M3A_PREREG.md)、結果は [VALIDATION.md](VALIDATION.md) の「M3a」にある。

## E7. 法線出射の単一損失スペクトルと有限の損失窓の面積（外から与える一定の長さ）

`simonsen1997_single_loss_spectrum_external_v0`（開発版、M4a）は、原著の式(19)のゼロ損失線と単一損失の衛星を、検出の運動エネルギーのビン積分として返す。有限の損失窓の面積と、窓の下限への感度も返す。核は E6 の深さ平均をそのまま使い、物理は足さない。[Simonsen ら (1997) 式(19)、p. 1615](https://doi.org/10.1103/PhysRevB.56.1612)

原著の式(19)は \(J(E)\propto\cos\theta\,[F(E)+\lambda\int_E^\infty dE'\,F(E')K_{\rm sc}(E_0,E'-E,\theta)]\) である。本実装は \(\theta=0\)、比例定数 1 とし、損失 \(T=E'-E\) の積分を有限の窓 \(W=[a,b]\) に限る。

\[
J_0(E)=F(E),\qquad
J_{1,j}(E)=L\int_W F(E+T)\,K_{{\rm sc},j}(E_0,T;L)\,dT,\qquad
J(E)=J_0(E)+J_{1,\rm total}(E)
\]

\[
A_j=\int_a^bK_{{\rm sc},j}(E_0,T;L)\,dT,\qquad
w_{1,W}=L\,A_{\rm total},\qquad
H_j=L\,a\,K_{{\rm sc},j}(E_0,a;L)
\]

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| \(F\)、\(I\) | 一次線と強度の単位。δ 線（位置 \(E_i\)、強度 \(s_i\)）か、ビン（端 \(e_m\)、積分強度 \(f_m\)）のどちらか | \(I\) |
| \(E\)、\(E'\)、\(T\) | 検出の運動エネルギー、一次励起の運動エネルギー、損失。\(E=E'-T\) | eV |
| \(E_0\) | 核を評価するエネルギー（設定の `electron_energy_eV`）。すべての \(E'\) で同じ値に固定する | eV |
| \(L\)（`weight_length_A`） | 外から与える一定の長さ。既定の規約では、式(18)の重みと式(19)の乗数に同じ値を使う。既定値はない。重みの長さと乗数を分ける規約は下の E7a | Å |
| \(W=[a,b]\)（`loss_window_eV`） | 損失窓。必須で、既定値はない。\(0<a<b<E_0\) | eV |
| \(K_{{\rm sc},j}\) | E6 の深さ平均の出力（8項、`no_hole`、`hole`、`total`） | Å⁻¹ eV⁻¹ |
| \(A_j\)（`window_area`） | 窓の面積（11列） | Å⁻¹ |
| \(w_{1,W}\)（`single_loss_weight`） | 単一損失の面積の比 | 無次元 |
| \(H_j\)（`cutoff_sensitivity`） | 窓の下限への感度（11列） | 無次元 |
| \(D\) | 出力の領域。出力の端 \(u_0<\dots<u_N\) のビン \([u_n,u_{n+1})\)（最後だけ右端を含む）と、両側の外 \((-\infty,u_0)\)、\((u_N,\infty)\) | eV |

- **規約:** \(J_0=F\) に足す衛星は `total` の列である。重みの \(L\) と乗数の \(L\) が同じ値になるのは、均質・半無限、一様な生成、法線出射、弾性散乱なし、直線の軌道と一定の速度、ほかに損失しない確率 \(e^{-x/L}\) という仮定の導出の中だけである（\(\int_0^\infty e^{-x/L}dx=L\)、\(\int_0^\infty xe^{-x/L}K_{{\rm eff},j}\,dx=L^2K_{{\rm sc},j}\)）。核は \(E_0\) に固定し、\(E'\) ごとに変えない。導出は [M4a の設計](M4A_DESIGN.md) §2.2。
- **窓が必須であること:** 原著は、式(19)の積分が発散するので Al の計算を \(E+0.05\) eV から積分し、下限への感度は小さかったと書く（p. 1618）。本実装は下限を既定値にせず、窓 \([a,b]\) を必ず呼び出す側に指定させる。窓の外の損失は評価せず、収支にも入れない。活動線のある材料で、窓がギャップ端をまたぐ場合とギャップ下の場合は拒否する（`window_not_above_gap`）。
- **\(w_{1,W}\):** 一次線の全強度 \(S=\int F\) が正のとき、全領域の衛星の和と \(S\) の比は \(w_{1,W}=LA_{\rm total}\) になる。確率ではなく、1 を超えうる。例えば ELF データ集の Al（`benchmarks/inputs/al_elfdb_v1.0.bulk.json`）、\(E_0=1000\) eV、\(L=20\) Å、窓 \([1,50]\) eV、`k_total` では \(A_{\rm total}=0.06273\) Å⁻¹、\(w_{1,W}=1.2545\) になる（例示。登録した比較の条件の外で、独立な照合はしていない）。
- **ビン積分:** 領域 \(D\) の値は \(J_0(D)=G_D(0)\)、\(J_{1,j}(D)=L\int_a^bK_{{\rm sc},j}(T)\,G_D(T)\,dT\) で、\(G_D(T)=\int\mathbf 1_D(E'-T)\,F(dE')\) は一次線のうち損失 \(T\) で領域 \(D\) に入る強度である。スペクトルに返すのは `total`、`no_hole`、`hole` の3成分で、8項は面積だけを返す。負の値は切らない。両側の外も同じ定義で直接に積分し、全体からの残差では作らない。密度（`*_per_eV`）はビンの値をビン幅で割ったものである。
- **収支:** \(\sum_DJ_0=S\)、\(\sum_DJ_{1,j}=SLA_j\)、\(\sum_DJ=S(1+LA_{\rm total})\)。面積でも各領域でも `total = no_hole + hole_charge * hole`。差は `diagnostics.balance` に残し、補正しない。
- **片側の感度:** 領域ごとの感度は \(H^\pm_{j,D}=L\,a\,K_{{\rm sc},j}(a\pm0)\,G_D(a\pm0)\)（3成分、`cutoff_sensitivity_bins`）。δ 線では \(G_D\) が \(T=a\) で跳びうるので、左右の極限を別々に返す。
- **逆数の区間:** `inverse_window_area_A` は \(1/A_{\rm total}\)、`inverse_window_area_interval_A` は区間 \([1/(A+e_A),\,1/(A-e_A)]\) である（\(e_A\) は面積の推定誤差。\(A\) と \(e_A\) を厳密な有理数として求め、下端は下へ、上端は上へ丸める）。\(A-e_A\le0\) なら `None` と理由 `nonpositive_or_unresolved_area`、表せなければ `inverse_not_representable` を返す（`inverse_reason`）。この量を平均自由行程とは呼ばない。
- **単位:** 面積は Å⁻¹、\(w_{1,W}\) と \(H_j\) は無次元、ビン・両側の外・\(H^\pm\) は一次線の強度の単位 \(I\)、密度は \(I\)/eV、逆数は Å。軸は運動エネルギーである。
- **数値の方法:** 損失の方向を、Chebyshev–Lobatto の 17 点の補間のパネルに分ける。パネルの幅は、強度が正の振動子の最小の線幅、低損失と上端の比 2 の刻み、計算の前に求めた候補の特異点（体積の極と、動く運動量の端の根）からの距離で決める。各パネルで別の 16 点の検証の点を評価する。出力ごとの推定誤差は、補間、検証の残差、核の点の誤差、端、丸めの5項の和である。方法は [M4a の設計](M4A_DESIGN.md) §3。
- **成功の意味:** 成功は「判定対象の全出力（面積 11 列と、各ビン・両側の外の \(J_0\) と衛星の3成分）で、推定誤差が予算の内」という意味である。予算は、面積が \(\max(\varepsilon_{{\rm abs},A},\varepsilon_{\rm rel}|A_j|)\)、スペクトルが \(\max(\eta S_{\rm abs},\varepsilon_{\rm rel}|y|)\)（\(S_{\rm abs}\) は強度の絶対値の和）。既定は \(\varepsilon_{\rm rel}=10^{-6}\)、\(\varepsilon_{{\rm abs},A}=10^{-10}\) Å⁻¹、\(\eta=10^{-9}\)、`panel_limit` 400、`max_cells` \(2\times10^6\)。E5・E6 と同じ基準で、証明された上界ではない。解像の尺度の根拠の強さは `diagnostics.resolution_basis` に記録する。`fixed_poles_enumerated`（全活動線が無分散でギャップなし、全根を確かめた）、`endpoint_roots_enumerated`（分散ありでギャップなし、端の全根を確かめた）、`computed_candidates_only`（ギャップのある材料、または確かめられない根や列挙の失敗があった実行）の3つである。
- **失敗:** 入力の誤りは `SingleLossInputError`（`reason` は `invalid_window`、`window_not_above_gap`、`invalid_primary`、`invalid_output_edges`、`invalid_settings`、`invalid_schema`、`cell_limit_exceeded`、`resolution_exceeds_panel_limit`、`panel_width_unresolved`）、未収束は `SingleLossConvergenceError`（`kernel_not_converged`、`kernel_error_exceeds_budget`、`loss_budget_not_met`、`cell_limit_reached`、`panel_width_unresolved`、`nonfinite_output`）である。必要なパネル数が `panel_limit` を超えるときは、核を評価する前に拒否する。
- **警告:** 常に出すのは、窓の外の損失を含まないこと、\(w_{1,W}\) は確率でないこと、核を \(E_0\) に固定すること、一次線は内因性の損失の衛星を含まないと仮定すること（スペクトルのとき）と、E6 から継承した警告である。条件つきで出すのは、\(b/E_0>0.1\)、\(a\) が一次線の最小のビン幅より小さい、\(E_0\) が正の強度の点の最小と最大の間にない、一次線の台の下端が \(b\) より小さい、`resolution_basis` が `computed_candidates_only`、照合した条件の外、正孔電荷が 1 でない、の各場合である。
- **照合の表示:** `run.json` の `comparison` は、配布する完了した比較の表 `src/xpsloss/data/m4a_completed_comparisons.json` と入力の条件（正規化した数値の列）を比べる。全条件が登録の記録と一致し、その群の必須キーがすべて合格のときだけ `matched`、核の条件だけが合えば `kernel_conditions_checked`、ほかは `not_matched` である。`physical_validity` は常に `not_established`。

**扱わないもの:** 多重散乱、斜出射、軌道ごとの長さと自己無撞着な長さ、弾性散乱、窓の外の損失、内因性の損失の分離、\(E'\) ごとに変わる核、装置関数と寿命幅の付加、背景の差引き、フィット。QUEELS-XPS の完全なモデルを実装したものではない。

| 式 | 実装 | 検証 |
| --- | --- | --- |
| \(A_j\)、\(w_{1,W}\)、\(H_j\)、逆数の区間 | `xps_single_loss.loss_window_area`（パネルは `_loss_panels`、候補の根は `_loss_roots`）。CLI は `xps-loss-window-area` | 解析の参照 R-A、実材料の参照 R-B と CE（[M4a の事前登録](M4A_PREREG.md)、[VALIDATION.md](VALIDATION.md) の「M4a」） |
| \(J_0\)、\(J_{1,j}\)、\(J\)、両側の外、\(H^\pm_{j,D}\) | `xps_single_loss.single_loss_spectrum`（セルの積分は `_fold`）。CLI は `xps-single-loss` | 同上 |

### E7a. 重みの長さと式(19)の乗数を分ける規約（M4a-2）

E7 の \(L\) を、式(18)の重みに使う長さ \(L_w\) と、式(19)の乗数 \(m\) に分ける。**物理は足さない。** 核（E6）、窓の積分、セルの積分、誤差の予算は E7 のままである。どの規約が物理として正しいかは主張しない。比較のために規約を選べるようにしたものである。

\[
J_{1,j}(E)=m\int_W F(E+T)\,K_{{\rm sc},j}(E_0,T;L_w)\,dT,\qquad
w_{1,W}=m\,A_{\rm total}(L_w),\qquad
H_j=m\,a\,K_{{\rm sc},j}(E_0,a;L_w)
\]

面積 \(A_j(L_w)\) と \(J_0\) は乗数によらない。収支は \(\sum_DJ_{1,j}=S\,m\,A_j\)、上端の診断 `LbK_b` は \(m\,b\,K_j(b)\) である。

| 重みの長さ（`weight_length_kind`） | \(L_w\) | 根拠 |
| --- | --- | --- |
| `external_constant`（既定） | 呼び出す側が与える正の値 `weight_length_A` | 式(18) [Simonsen ら (1997) p. 1615](https://doi.org/10.1103/PhysRevB.56.1612) |
| `inverse_13a_window_area` | \(1/D\)、\(D=A_{13a}(W)+c\)。`weight_length_A=None` を明示する | 文献の式ではない。面積と \(\lambda_c^{-1}\) の和の逆数という形（[Tougaard ら (2022) p. 4 式(7)](https://doi.org/10.1002/sia.7095)。面積は 0 から \(E_{\max}\) までの全有効断面積）を、13a の項と有限の窓に当てはめたもの。根拠は M3a-2 の判定（[ROADMAP](ROADMAP.md)） |

| 乗数（`multiplier`） | \(m\) | \(w_{1,W}\) | 根拠 |
| --- | --- | --- | --- |
| `weight_length`（既定） | \(L_w\) | \(L_wA_{\rm total}\) | 式(18)–(19)（同じ \(\lambda\)） |
| `external` | 呼び出す側が与える正の値 `multiplier_A`（Å） | \(mA_{\rm total}\) | 文献の式ではない（比較のための組み合わせ） |
| `inverse_window_area`（S2） | \(1/A_{\rm total}(W;L_w)\) | 定義により 1 | 深さ平均の後の断面積の面積で乗数を決める規約（[Yubero と Tougaard (2005) p. 045414-3 式(9)の後](https://doi.org/10.1103/PhysRevB.71.045414)）を、有限の窓に限って読んだもの。文献の面積は計算した損失の全範囲である |

- **規約の組と識別子:** 外部の重みで乗数が \(L_w\)（省略、`weight_length`、または \(L_w\) と binary64 で等しい `multiplier_A`）のときは、E7 と同じ演算の列でビット単位で同じ値を返し、識別子 `simonsen1997_single_loss_spectrum_external_v0`、照合、警告も E7 のままである。それ以外の5通りの組は識別子 `simonsen1997_single_loss_spectrum_conventions_v0`（式(19)の形に別の規約の長さを入れたもの）で、照合は常に `not_matched` である。組の全体は [M4a-2 の設計](M4A2_DESIGN.md) §2.4。
- **\(L_w\neq m\) は導出の外:** E7 の導出で重みの長さと乗数が同じになるのは、一様な生成と一定長の指数の減衰を仮定したときだけである。\(L_w\neq m\) の組はその導出から外れる規約の組み合わせで、導出の結果ではない（警告 W-len）。
- **窓の 13a の長さ:** \(c\ge0\)（`extra_inverse_length_per_A`、Å⁻¹、既定 0）は呼び出す側が与える追加の減衰で、パッケージは \(\lambda_c\) を計算しない。13a の被積分関数に深さも長さも入らないので、\(A_{13a}\) は重みの長さによらない。この長さは**平均自由行程ではない**。窓に依存し（ELF データ集の Al、1000 eV、`k_total` で、上端 50 eV と 99.8 eV でおよそ 27.06 Å と 25.96 Å）、軌道ごとの重みに近いことを確かめたのは、M3a-2 の列挙した条件（Al と Si、1000 eV、法線出射、正孔電荷 1）だけである（警告 W-13a、\(c>0\) なら W-c）。`k_parallel` の 13a は無限媒質のバルクの値と違うので、名前に bulk を使わない。外部から文献の IMFP を入れる使い方は、外部の重みで行う。
- **2段の計算:** 窓の 13a の重みでは、まず種の長さ（`weight_length_seed_A`、省略で 20 Å）で面積を積分し、`term_13a` だけを判定する（段α）。\(A^\alpha_{13a}\)、その誤差 \(e^\alpha\)、\(c\) を厳密な有理数で足し、\(L_w\) を \(1/D\) の最も近い binary64、区間 `weight_length_interval_A` を \([1/(D+e^\alpha),1/(D-e^\alpha)]\) の外向きの丸めとする。次に \(L_w\) で本計算をする（段β）。最後に \(|A^\alpha_{13a}-A^\beta_{13a}|\le e^\alpha+e^\beta\) を有理数で確かめ、外れたら `weight_length_inconsistent`。例外は `phase`（`weight_length_seed`、`production`）を持つ。既定の組の例外、計算の前の検査の例外、段αを飛ばした構造的なゼロの例外では `phase` は `None` である。
- **分母の規則:** 逆数を取る量 \(D\)（段αは \(A^\alpha_{13a}+c\)、S2 は \(A_{\rm total}\)）には、E7 の面積の予算に加えて \(D-e>0\) と \(e\le\varepsilon_{\rm rel}|D|\) を課す。満たさない段では逆数を保留し、畳み込みと成功をしない。保留の間は、分母の列の局所の誤差が最大のパネルを**1枚だけ**二分する（同じなら長いもの、さらに同じなら左のもの。通常の分割の候補は使わない）。符号は条件を満たした段で初めて決め、負なら `weight_length_undefined`、`multiplier_undefined`。分割が `panel_limit` で実行できないときに分母が決まっていなければ `weight_length_unresolved`、`multiplier_unresolved`。逆数か区間の端が表せなければ `*_not_representable`。強度が正の振動子がない材料は、核を呼ぶ前に判別する（窓の 13a の重みで \(c=0\) と、外部の重みで S2 は拒否。\(c>0\) なら \(L_w=1/c\) で段αを飛ばす）。
- **S2:** 段βの各段で面積を求めた後、分母が決まった段で \(m\) を決め直し、全パネルをその \(m\) で畳み直す。\(w_{1,W}\equiv1\) なので衛星の相対面積は予測ではなく定義で、意味をもつのは窓の中の形とビンごとの配分だけである（警告 W-S2）。一次線の強度が 0 のとき、\(w=1\) は規約の値で、衛星と一次線の比は定義されない。ポアソンの数え方による乗数は文献の式ではないので実装しない（一定長では \(L_w\) と一致する）。
- **感度:** `cutoff_sensitivity`、`cutoff_sensitivity_bins`、`LbK_b` は、報告した \(L_w\) と \(m\) を**固定した**偏感度である。窓を変えたときの \(m\) や \(L_w\) の変化は含まない（S2 では \(w\equiv1\) なので全微分は 0 だが、\(H\) は 0 でない。警告 W-sens）。窓を変えたときの値が要るなら、その窓で計算し直す。
- **返す誤差:** 乗数に区間があるとき（S2 と、窓の 13a の重みで乗数 `weight_length`。`multiplier_interval_A`）、乗数に比例する出力の誤差を \(e_{\rm ret}=e_f+\rho(|\hat y|+e_f)\)、\(\rho=\max(\hat m-m_-,m_+-\hat m)/\hat m\) とする（厳密な有理数で求め、最後に1回だけ上向きに丸める。`diagnostics["multiplier"]["rho"]`）。成功の判定と `diagnostics.error_budget` は、報告した \(L_w\) と \(\hat m\) を固定した \(e_f\) で行う。そのため、この場合の `*_error` は `error_budget` の `estimate` と一致せず、\(\varepsilon_{\rm rel}\) が小さいとき予算のおよそ2倍までになりうる。S2 の \(w\) の誤差は \(|w-1|\) の上向きの丸めである。**重みの長さの不確かさが核 \(K(L_w)\) を通して出力を動かす分は、誤差に入れない。** その大きさは区間 `weight_length_interval_A` で示す（M3a-2 の条件では、長さ ±10 % で面積が数 %）。
- **検査の順:** 新しい引数（\(c\)、種、乗数）の検査は、E7 の窓・材料・設定の検査の後、出力の端と一次線の検査の前に行う。どれも核を評価する前の拒否である。

| 式 | 実装 | 検証 |
| --- | --- | --- |
| \(L_w\)、\(m\)、区間、2段の計算、分母の規則、返す誤差 | `xps_single_loss.loss_window_area`、`single_loss_spectrum`（キーワード `weight_length_kind`、`extra_inverse_length_per_A`、`weight_length_seed_A`、`multiplier`、`multiplier_A`）。分母と再畳み込みは `_loss_panels`、返す誤差は `_fold`。CLI は `--weight-length-inverse-13a-window-area`、`--extra-inverse-length-per-A`、`--weight-length-seed-A`、`--multiplier-A`、`--multiplier-inverse-window-area` | R-B から写像した参照（[M4a-2 の事前登録](M4A2_PREREG.md)、[VALIDATION.md](VALIDATION.md) の「M4a-2」） |
