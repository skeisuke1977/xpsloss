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

## E8. 法線出射の多重損失（主模型 `first_loss_13a_chain_v0`）

`first_loss_13a_chain_v0`（開発版 0.4.0.dev0、M4b-1）は、E7・E7a のゼロ損失線と単一損失に、同じ損失窓の中の損失を繰り返した次数 \(n\ge2\) の寄与を足す。**近似の模型である。** 最初の損失には E7a の全核（11 列）を使い、2回目以降の損失には 13a の項だけを使う。核（E6）は足さない。軸は運動エネルギーだけである（束縛エネルギーの入出力は M4b-2）。設計は [M4b の設計](M4B_DESIGN.md)（第5版）、検証の条件は [M4b-1 の事前登録](M4B1_PREREG.md) にある。

\[
K_{1,j}(T)=m\,\mathbf 1_W(T)\,K_{{\rm sc},j}(E_0,T;L_w),\qquad
R(T)=m\,\mathbf 1_W(T)\,K_{13a}(E_0,T)
\]

\[
C_0=\delta,\qquad C_{1,j}=K_{1,j},\qquad C_{n,j}=K_{1,j}*R^{*(n-1)}\ (n\ge2)
\]

\[
J_{n,j}(E)=\int C_{n,j}(s)\,F(E+s)\,ds,\qquad
J^{(N)}(E)=J_0(E)+\sum_{n=1}^{N}J_{n,\rm total}(E)
\]

\[
w_{n,j}=m\,A_j\,r^{\,n-1},\qquad r=m\,A_{13a}
\]

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| \(W=[a,b]\)、\(L_w\)、\(m\)、\(F\)、\(E_0\)、\(A_j\) | E7・E7a と同じ（損失窓、重みの長さ、乗数、一次線、核を評価するエネルギー、窓の面積） | — |
| \(K_{1,j}\) | 1回目の損失の核。\(j\) は 11 列（8項、`no_hole`、`hole`、`total`） | eV⁻¹ |
| \(R\) | 繰り返しの核。13a の列だけ（`repeated_kernel` は `term_13a`） | eV⁻¹ |
| \(C_{n,j}\) | 次数 \(n\) の損失の分布。台は \([na,nb]\)。\(s\) は累積の損失 | eV⁻¹ |
| \(J_{n,j}\) | 次数 \(n\) の寄与。領域ごとのビン積分で返す（`orders[n-1]` の `total`、`no_hole`、`hole`） | \(I\) |
| \(J^{(N)}\)（`spectrum`） | 次数 \(N\) までの和 | \(I\) |
| \(w_{n,j}\)（`area`） | 次数 \(n\) の全軸の面積の比（11 列）。確率ではない | 無次元 |
| \(r\)（`order_area_ratio`） | 隣り合う次数の面積の比 | 無次元 |
| \(N\)（`max_order`） | 計算する最大の次数。必須で、既定値はない | — |

- **根拠:** 原著は、式(18)の核が最初の非弾性散乱だけを表し、分けられない内因性の部分を含むので、繰り返しの畳み込みで全スペクトルを作れないと書く [Simonsen ら (1997) p. 1615](https://doi.org/10.1103/PhysRevB.56.1612)。そこで全核は1回目にだけ使う。全核をそのまま繰り返す文献の式は、各次数の係数が文献の間で食い違うので、公開の模型にしない（[M4b の設計](M4B_DESIGN.md) §2。式(19)と多重散乱の式の原典の Tougaard と Sigmund (1982) などは論文庫になく、その記述は二次文献による）。
- **正孔・干渉・表面・真空は高々1回:** 2回目以降に入るのは 13a だけである。13b、14、15-T1、15-T2、15-T3、16-e、16-h は1回目の核にだけ入る（`hole_treatment` は `first_event_only`）。\(\bar g=K_{1,\rm total}-R\) とおくと \(\sum_nC_{n,\rm total}=(\delta+\bar g)*\sum_{k\ge0}R^{*k}\) で、\(\bar g\) は軌道に高々1回入る。繰り返しの 13a は、1回目と同じ運動量の指定の列である。分散のある `k_parallel` の 13a は無限媒質のバルクの値と一般に違うので、警告を出す。
- **各次数の係数:** 半無限・一様な生成・法線出射・直線の軌道・一定の速度・弾性散乱なしで、ほかに損失しない確率を \(e^{-x/L}\) とすると、次数 \(n\) の係数は \(L^n\) になる（[M4b の設計](M4B_DESIGN.md) §3）。そこで繰り返しの因子にも1回目と同じ乗数 \(m\) を掛け、次数 \(n\) は \(m^n\) に比例する。\(L_w\) は1回目の核の深さの重みにだけ入る。\(m\ne L_w\) の組はこの導出の外にある（E7a の W-len を引き継ぐ）。
- **核は \(E_0\) に固定:** 全次数で同じ \(E_0\) の核を使い、次数ごとにエネルギーを下げない（`kernel_energy_policy` は `fixed_E0`）。\(Nb/E_0>0.1\) なら警告する。
- **成分と符号:** \(C_n\) は \(K_1\) に線形なので、各次数・各領域・面積で `total = no_hole + hole_charge * hole` が成り立つ。非負の \(F\)、\(K_{1,\rm total}\)、\(R\) から作る `total` は全次数で非負である。`hole` は負になりうる。負の値は切らない。`total` と 13a の核の、評価点の最小と連続区間の最小を `diagnostics["negative_values"]` に残し、負なら警告する。
- **面積:** 次数 \(n\ge2\) の \(w_{n,j}\) は、その次数のパネルの面積から求める。恒等値 \(w_{1,j}\hat r^{\,n-1}\) で上書きしない。`order_area_ratio` は \(\hat r=\hat mA_{13a}\) とその誤差・区間である。`convention_ratio` は、規約から決まる \(r\) の性質を別に示す（窓の 13a の重みと乗数 `weight_length` の組で \(c=0\) なら `exactly_one`、\(c>0\) なら `below_one`。ほかは `not_determined_by_convention`）。`area_zero_kind` は、面積が 0 の列を `exact_zero`（構造の証拠がある）と `rounded_zero`（中心の値が 0 に丸まっただけ）に分ける。今の実装が `exact_zero` とするのは、強度がゼロの材料の場合だけである。
- **窓の下限への感度:** \(h_{1,j}=m\,a\,K_{{\rm sc},j}(a)\)（E7・E7a の \(H_j\)）、\(n\ge2\) で \(h_{n,j}=h_{1,j}\,r^{\,n-1}+(n-1)\,w_{1,j}\,h_r\,r^{\,n-2}\)、\(h_r=h_{1,13a}\)。\(m\)、\(L_w\)、\(b\)、\(E_0\)、正孔電荷を固定した偏感度である（E7a の W-sens）。ビンごとの片側の感度は1次だけ（E7 のまま）。
- **`max_order` は必須:** 整数 \(N\ge1\) か `"display_support"` を与える。既定値はない。正の強度の一次線の台 \([p_{\min},p_{\max}]\) と出力の下端 \(d_{\min}\) から、表示のビンに届く次数の上限は \(N_D=\max(0,\lfloor(p_{\max}-d_{\min})/a\rfloor)\)（`truncation.support_order`。入力の binary64 を厳密な有理数として計算し、等号の次数は残す）。`"display_support"` は \(N=N_D\) とする。面積だけの入口は整数だけを受け付ける。
- **打ち切りの状態（`truncation.state`）:** \(N\ge N_D\) なら `complete_on_display`。そうでなく、`neglect_tolerance_fraction`（既定値なし）を与え、表示のビンすべての上界が許容 \(\times S_{\rm abs}\) 以下なら `bounded_on_display`。許容を与えなければ `truncated` で、成功として返す（落とした次数は上界つきで返すだけで、警告を出す）。許容を与えて上界が超えれば `truncation_bound_exceeds_tolerance`（必要な次数の見積りつき）、上界が求まらなければ `truncation_bound_unavailable` で止まる。両側の外は判定に入れない。
- **落とした次数の上界（領域ごとの累積）:** 領域 \(D\) の下端を \(d_D\)、\(U_D=p_{\max}-d_D\)、\(\Lambda_f(\tau)=\int_We^{-\tau t}|f(t)|\,dt\) とすると、\(\Lambda_R(\tau)<1\) の \(\tau>0\) で、次数 \(>N\) のうち損失 \(\le U_D\) に落ちる絶対の量は

\[
S_{\rm abs}\,e^{\tau U_D}\,\Lambda_{K_1}(\tau)\,\frac{\Lambda_R(\tau)^N}{1-\Lambda_R(\tau)}
\]

  以下である。`truncation.neglected_at_or_above_region_bound` は、この式の \(\tau\) の格子（\(\tau_i=10^{-3+i/20}\) eV⁻¹ のうち \(\tau_i\le40/a\) のもの。\(\|R\|_1\) の上側の値が 1 未満なら \(\tau\to0\) も候補）での最小である。**ビンごとの値ではなく、「領域 \(D\) とそれより高い運動エネルギーの側の全部に落ちる量」の累積の上界である。** \(\Lambda_f\) はパネルごとの上側の値（乗数の区間の上端と、パネルの誤差を含む）で作り、対数で計算して外向きに丸める。上流が推定の誤差なので、推定の誤差を入れた上界であり、数学的な保証ではない。\(N\) より高い次数が届かないビンは 0（`bound_reason` は `complete`）、使える \(\tau\) がなければ `None`（`no_usable_tau`）、表せなければ `None`（`not_representable`）。
- **落とした次数の全軸の面積（`truncation.neglected_total_area`）:** \(S=0\) か1回目の核がゼロなら 0。\(R=0\) なら \(N\ge1\) で 0（\(N=0\) では \(Sw\)）。繰り返しの核が非負なら、\(\hat r\) の区間の上端が 1 未満のとき \(S\,w\,r^N/(1-r)\) とその区間、下端が 1 以上なら `None`（`diverges`）、区間が 1 をまたげば `None`（`undetermined`）。繰り返しの核に負の値があれば `None`（`repeated_kernel_not_nonnegative`）。両側の外の値には、落とした次数を足さない。
- **数値の方法:** 1次は E7・E7a の計算を同じ引数で1回だけ行い、その最終のパネル（17 点の補間）を受け取る。次数 0 と 1 の出力は E7・E7a の結果そのもので、ビット単位で同じである（結果の `single_loss`）。次数 \(n\ge2\) は \(C_n(s)=\int C_{n-1}(s-t)R(t)\,dt\) を、両方のパネルの境で区切った区間ごとに 17 点のガウス・ルジャンドル則で積分し、17 点の補間のパネルに入れ直す。外側の網は、窓の端の \(n\) 個の和 \(ka+(n-k)b\)（\(k=0..n\)）を必須の節にして始め（有理数で作り、重複を除いて1回だけ丸める）、二分で細かくする。領域への畳み込みは E7 と同じ幾何とセルの規則で、各次数の最終の網に行う。乗数は有理数の \(\hat m^n\) を最後に1回掛ける。
- **誤差の3つの分け方:** 次数 \(n\) の誤差を、(i) 1次のパネルから継承した分（`inherited`）、(ii) \(2..n-1\) 次の網で生じて伝わった分（`propagated`）、(iii) \(n\) 次の網で生じた分（`local`）に分けて返す（`*_error_parts`）。\(n\) 次の網の二分で減るのは (iii) だけである。
  - (iii) は `order_epsrel`（\(\varepsilon_o\)）、`order_epsabs_fraction`（\(\eta_o\)）で判定する。予算は、ビンが \(\max(\eta_oS_{\rm abs},\varepsilon_o|y|)\)、面積が \(\max(\eta_o,\varepsilon_o|w_n|)\)。満たすまで二分し、上限に当たれば `order_budget_not_met`。
  - (i) は1次の設定 `epsrel`（\(\varepsilon\)）、`spectrum_epsabs_fraction`（\(\eta\)）、`area_epsabs`（\(\eta_A\)）で判定する。\(\hat w_{{\rm abs},j}=\|\hat K_{1,j}\|_1\)、\(\hat r_{\rm abs}=\|\hat R\|_1\)、\(\phi_{n,j}=\hat r_{\rm abs}^{\,n-1}+(n-1)\hat w_{{\rm abs},j}\hat r_{\rm abs}^{\,n-2}\) として、予算は、ビンが \(\max(\phi_{n,j}\,\eta\,S_{\rm abs},\ n\varepsilon|y|)\)、面積が \(\max(\phi_{n,j}\,\hat m\,\eta_A,\ n\varepsilon|w_{n,j}|)\)。1次のパネルは作り直さないので、超えたら `inherited_error_exceeds_order_budget` で止まる（次数、超えたキー、相対と絶対のどちらの側か、下げる設定の名前つき）。
  - (ii) は実額を返す。合否の判定は課さない。
  - **既定の `epsrel`（\(10^{-6}\)、E7 と同じ）では、`max_order` が 2 以上の典型の実行（Al、δ の一次線、1 eV のビン。滑らかな解析の核でも）が、2次で `inherited_error_exceeds_order_budget` で止まる。`epsrel=1e-7` を与えると通る。** 既定値は E7 と同じに保ち、例外の文で `epsrel=1e-7` を案内する（利用者の決定、2026-11-27）。
  - 高次の誤差の内訳と、乗数を含む集計は、有理数で足して外向きに丸める。パネルに由来する誤差そのものは E7 と同じ推定で、証明された上界ではない。各次数の判定を通っても、伝播の分と乗数の区間で、高次の返す誤差は1次の設定より大きくなりうる。
- **乗数の区間:** 乗数に区間 \([m_-,m_+]\) があるとき（E7a）、\(\rho_n=\max(1-(m_-/\hat m)^n,\ (m_+/\hat m)^n-1)\) として、返す誤差を \(e_f+\rho_n(|\hat y|+e_f)\) とする（\(e_f\) は固定の乗数 \(\hat m\) での (i)–(iii) の和。上へ丸める）。再帰には固定の乗数の誤差だけを渡し、区間は各次数の最終の出力に1回だけ当てる（内訳の `multiplier_interval`）。\(L_w\) の区間が核を通して出力を動かす分は、E7a と同じく誤差に入れない。
- **和の誤差:** `spectrum_error` は \(J_0\) と各次数の返す誤差の絶対和に、加算の丸めを足したものである。和そのものに予算は課さない。
- **資源の上限:** `order_limit`（既定 64）、`max_convolution_points`（既定 \(2\times10^8\)。積分した区間の数 \(\times17\) の、2次から \(N\) 次までの累計）、`max_cells`（既定 \(2\times10^6\)、次数ごと）。入力だけで決まる超過は核を評価する前に `MultipleLossInputError`（`order_limit_exceeded`、`work_lower_bound_exceeded`）、1次の後の見積りの超過は `convolution_work_limit_exceeded`、適応の途中の超過は `convolution_work_limit_reached`。資源の上限であって、物理の既定値ではない。
- **失敗:** 入力の誤りは `MultipleLossInputError`（`reason` は `invalid_model`、`invalid_max_order`、`invalid_neglect_tolerance`、`invalid_settings`、`order_limit_exceeded`、`work_lower_bound_exceeded`、CLI の `invalid_schema`）、未収束は `MultipleLossConvergenceError`（`order_budget_not_met`、`inherited_error_exceeds_order_budget`、`order_mesh_unresolved`、`convolution_work_limit_exceeded`、`convolution_work_limit_reached`、`truncation_bound_exceeds_tolerance`、`truncation_bound_unavailable`、`nonfinite_output`。`order`、`stage`、`inner` を持つ）。E7・E7a の例外はそのまま通す。E7 の入力の検査（材料、運動量、角度、エネルギー、長さと窓、数値の設定、規約、出力の端）を同じ順に行い、その後、核を評価する前に、模型・次数・資源を調べる。
- **構造的なゼロ:** 一次線の強度がすべて 0 なら、ビンと両側の外は厳密に 0 で、面積は計算する。強度がゼロの材料では、高次の値・誤差・面積は厳密に 0 で \(J=J_0\)（E7a の未定義の規約の拒否が先）。
- **警告:** 常に出すのは、模型が近似であること、表面・正孔・干渉は1回目にだけ入ること、軌道ごとの深さの相関を落とすこと、弾性散乱を含まないこと、次数の面積は確率でないことと、E7・E7a から継承した警告である。条件つきで出すのは、状態が `truncated`、\(\hat r\) の区間の上端が 1 以上、\(m\ne L_w\)、\(Nb/E_0>0.1\)、`total` か 13a の核に負の値、分散のある `k_parallel`、\(E_0>5\) keV、の各場合である。
- **照合の表示:** `run.json` の `comparison` は、群 `orders`、`areas`、`spectrum` が常に `not_matched` である（完了した比較の表は M4b-1 では作らない）。`physical_validity` は常に `not_established`。数値の検証が確かめたのは、次数ごとの畳み込みが、与えた核と模型の定義に対して正しいことだけで、模型の物理の妥当性は確かめていない。
- **近似の大きさ（模型の不確かさ、設計の試算）:** 表面・正孔などの超過もポアソンで数える模型との差は、2次の全軸の面積で Al 15.0 %、Si 17.4 %。軌道ごとの深さの相関を切ることによる差は、交差項の面積では 1 % 未満だが、形では同じ桁になりうる（[M4b の設計](M4B_DESIGN.md) §4.4。設計の判断のための試算で、検証した値ではない）。

**扱わないもの:** 束縛エネルギーの入出力（M4b-2）、実測スペクトルとの比較（M4b-3）、斜出射、弾性散乱、軌道ごとの長さと自己無撞着な長さ、窓の外の損失、内因性の損失の分離、\(E'\) や次数ごとに変わる核、装置関数の付加、背景の差引き、フィット。全核の繰り返しなどの比較用の模型は非公開で、API・CLI・`capabilities` に出さない。QUEELS-XPS の完全なモデルを実装したものではない。

| 式 | 実装 | 検証 |
| --- | --- | --- |
| \(C_{n,j}\)、\(J_{n,j}\)、\(J^{(N)}\)、打ち切りの状態と上界 | `xps_multiple_loss.multiple_loss_spectrum`（高次の網・誤差の内訳・判定は `_loss_orders`、区間ごとの求積と必須の節は `_order_quadrature`、上界・感度・区間は `_order_bounds`、次数ごとの結果と和は `_order_results`、結果と例外の型は `_order_types`）。1次は `xps_single_loss` の計算を1回呼ぶ。CLI は `xps-multiple-loss`（`--model`、`--max-order` は必須） | 解析の参照 R-A、実材料の参照 R-B、期待する状態の 109 行、変異 29 件（[M4b-1 の事前登録](M4B1_PREREG.md)、[VALIDATION.md](VALIDATION.md) の「M4b-1」） |
| \(w_{n,j}\)、\(h_{n,j}\)、\(\hat r\) と区間 | `xps_multiple_loss.multiple_loss_order_areas`（一次線と出力の端を取らない。整数の `max_order` だけ） | 同上 |
