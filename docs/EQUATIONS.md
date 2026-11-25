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
