# M4b-1 の事前登録　検証の条件・許容幅・期待する結果

作成日: 2026-11-27。状態: **凍結（利用者の了解、2026-11-27。第4稿を凍結した）。**実装も参照値もない段階で固定した。Codex（実装担当）が第1稿（判定「修正後に凍結」。blocker 4、major 5、minor 1）、第2稿（同じ判定。blocker 2、major 5、minor 1）、第3稿（同じ判定。major 3）を確かめた（§10）。この第4稿は3回目の指摘を反映した（入力と許容幅は変えていない）。利用者が了解し（2026-11-27）、統括担当が凍結した。設計 §14 の D2 の扱い（継承の判定を超えたら止め、相対の側なら `epsrel` を下げるよう理由に書く。符号つきの成分だけ別の尺度で測らない）も、利用者が凍結の前に選んだ。この文書は、実装と参照値の計算の前に、検証の条件、許容幅、期待する結果を固定する契約になる。

**根拠:**

- [M4b の設計](M4B_DESIGN.md)（第5版）の §4–§6、§9、§10、§14。利用者が §12 の U1–U13 を推奨どおり決めた（2026-11-27）。
- 試走 T1–T13 の報告 `work/consult/m4b/trial/trial_report.md`。**試走の数値は期待値に使わない。**この文書が試走から取ったのは、既定値、具体的な入力、期待する終了状態、時間の見積りだけである。試作のコードは読んでいない。
- Codex の確認 `work/consult/m4b/prereg/codex_prereg_review.md`（1回目）、`codex_prereg_review2.md`（2回目）、`codex_prereg_review3.md`（3回目）。非公開の名前、結果の欄の形、注入の条件、変異の検出の経路は、この確認の案で決めた。採らなかった点と変えて採った点は `work/consult/m4b/prereg/open_questions.md` の「Codex の確認への返事」（1回目、2回目、3回目）に書いた。
- 下書きの担当が決めた細部は §10 の「ここで決めた事項」に並べた。まだ決まっていない点は同じ `open_questions.md` にある。

**範囲（U11 (a)）:** 畳み込みの数値、主模型 B1 `first_loss_13a_chain_v0`、運動エネルギー（KE）の API・CLI・`run.json`・`capabilities`、非公開の参照の模型。束縛エネルギーの入出力（設計 §7、§10 の「BE 軸」の群、変異 X25–X27・X32、理由 `invalid_energy_axis`・`axis_conversion_not_injective`・`axis_conversion_changes_support`）は M4b-2 に回す。実測との比較は M4b-3。法線出射だけ。QUEELS-XPS の完全モデルは対象でない。

**凍結の後の規則:**

- この文書と期待する状態の表を含むコミットの後に、参照値の計算と実装を始める。
- 比較の後に、許容幅と期待する結果を変えない。変える場合は、利用者の了解を得て、旧版を残し、別の試験として記録する（§11）。
- 必須の検査が時間の上限（§8）に収まらなければ、§5.3 の順で全件に回す。それでも収まらなければ旧版を残して再登録する。
- 不合格なら、(1) 1次（M4a の結果とのビット一致）、(2) 高次の畳み込み（点の値。層C と解析の入口）、(3) 網と二分、(4) 高次の畳み込みの後の領域への畳み込み（層B の高次版）、(5) 誤差の3分けと判定、(6) 打ち切りと上界、(7) 駆動・API・CLI、(8) 参照の経路、の順に切り分けて記録する。

記号は設計の「表記」と [M4a の事前登録](M4A_PREREG.md) に従う。\(K_j\) は 11 列（`term_13a`、`term_13b`、`term_14`、`term_15_T1`、`term_15_T2`、`term_15_T3`、`term_16_e`、`term_16_h`、`no_hole`、`hole`、`total`）。\(K_{1,j}=m\mathbf 1_WK_j\)、\(R=m\mathbf 1_WK_{13a}\)、\(C_{n,j}=K_{1,j}*R^{*(n-1)}\)、\(J_{n,j}=\mathcal A_{C_{n,j}}F\)、\(w_{n,j}=mA_jr^{n-1}\)、\(r=mA_{13a}\)。\(\varepsilon\)、\(\eta\)、\(\eta_A\) は1次の `epsrel`、`spectrum_epsabs_fraction`、`area_epsabs`、\(\varepsilon_o\)、\(\eta_o\) は `order_epsrel`、`order_epsabs_fraction`。\(\epsilon=2^{-52}\)。\(u_r\) は参照値の不確かさ。\(I\) は強度の単位。

## 1. 入力

数値は十進で与え、binary64 に変換して渡す。参照の計算にも同じ binary64 を渡す。機械で読める同じ内容（十進と16進）が `benchmarks/numerics/m4b1/expected_states_m4b1.json` の `materials`、`kernels`、`primaries`、`output_edges` にある。

### 1.1 材料

| 記号 | 内容 \((E,A,\Gamma,\alpha)\) | 来歴 |
| --- | --- | --- |
| Al | `benchmarks/inputs/al_elfdb_v1.0.json`（SHA-256 `17ded0465a5e42124eebf44d9573369a5cbea87df59b5df11216ea2f7432d1c4`）の材料。\(g=0\)、\((15.6,243.3598,1.5,1.0)\)。同じ数値の実行用の設定は `al_elfdb_v1.0.bulk.json`（`978931c1651881e6c7dfa739f4d626e4496225bf24bceb815466ff651934eae2`） | ELF データ集 v1.0（CC BY 4.0） |
| Si | `benchmarks/inputs/si_elfdb_v1.0.json`（`49020931fcb6d83f799b6a8b0cfd0a6eedb8c3a5dc446b88793bc9e4e5704079`）の材料。\(g=1.12\)、\((10,6.08458,5,0.5)\)、\((14,30.4229,5,0.5)\)、\((16.8,212.9603,3.8,0.5)\)。実行用は `si_elfdb_v1.0.bulk.json`（`471e513416369eec3dd12302a8e00823113935f78988293205b34f469f7200be`） | 同上 |
| AlZ | Al の強度を 0 にしたもの（構造的なゼロ） | synthetic |
| `synthetic_wide` | \(g=0\)、\((20,200,8,0)\)。CI 用の安い材料（M4a と同じ）。無分散 | synthetic |

表の生成は、記録のファイルと実行用の設定の材料の数値が等しいことを検査する。

### 1.2 共通の条件

- **公開の入口:** \(E_0=1000\) eV、外部の長さ \(L_w=20\) Å、\(Q_h=1\)、`k_total`、出射角 0°、乗数は `weight_length`（\(m=L_w\)）。窓は Al と `synthetic_wide` が \([1,50]\) eV、Si が \([5,50]\) eV（M4a の RB1–RB4 と同じ。1次の列を M4a の結果とそのまま比べられる）。
- **解析の入口:** \(L_w=2\) Å、\(Q_h=1\)、\(m=L_w\)。\(E_0\) は持たない。
- 行ごとに変えた点は §5.2 と表に書く。

### 1.3 解析の核（単位は Å⁻¹ eV⁻¹。書いていない列は 0。`no_hole`、`hole`、`total` は定義どおりに作る）

| 記号 | 窓（eV） | 定義 | 申告 |
| --- | --- | --- | --- |
| K_A1 | \([1.3,7.7]\) | 13a \(=0.05\)（箱形。\(K_1=R\)） | 跳び 1.3、7.7 |
| K_A1s | \([2,5]\) | 13a \(=0.1\)（端が二進で厳密） | 跳び 2、5 |
| K_A2 | \([1.3,7.7]\) | 13a \(=0.08e^{-0.25T}\)（傾けた箱） | 跳び 1.3、7.7、変動の尺度 4 eV |
| K_A3 | \([1,8]\) | 13a \(=0.04\)（\([2,5]\)）、15-T1 \(=0.03\)（\([1,3]\)）、14 \(=0.02\)（\([3,6]\)）、15-T2 \(=-0.05\)（\([1.5,2.5]\)） | 跳び 1.5、2、2.5、3、5、6 |
| K_A3z | \([1,8]\) | K_A3 の 15-T2 を \(-0.06\) にしたもの（`hole` の代数の面積が 0） | 同じ |
| K_A4_9、K_A4_16 | \([2,5]\) | 13a \(=0.05\,(1+P_d(u))\)、\(u=(2T-7)/3\)、\(d=9\)、16 | 跳び 2、5 |
| K_A5 | \([2,12]\) | 13a \(=0.05/[1+((T-5.3)/0.4)^2]\) | 特異点（5.3、0.4） |
| K_T1 | \([1,2.5]\) | 13a \(=0.2\)（\(r=0.6\)） | 跳び 1、2.5 |
| K_T2 | \([1,50]\) | 13a \(=0.0072\)、15-T1 \(=0.005\)（\(m=2\) で \(r=0.7056\)、\(w=1.1956\)） | 跳び 1、50 |
| K_T4_one、K_T4_above、K_T4_below | \([1,5]\) | 13a \(=1/8\)、\((1+2^{-10})/8\)、\((1-2^{-10})/8\)（\(m=2\) で \(r=1\)、\(1\pm2^{-10}\)。二進で厳密） | 跳び 1、5 |
| K_T5 | \([1,50]\) | 13a \(=0.01\) | 跳び 1、50 |
| K_SGN | \([2,5]\) | 13a \(=0.10\)、14 \(=+0.04\)（\([2,3.5]\)）、\(-0.04\)（\([3.5,5]\)） | 跳び 2、3.5、5 |
| K_LIN | \([2,5]\) | 13a \(=0.10\)、14 \(=0.02(T-3)\) | 跳び 2、5 |
| K_STEP | \([2,5]\) | 13a \(=k/16\)（\([2+0.5(k-1),\,2+0.5k]\)、\(k=1..6\)。1次は 6 枚、2次の必須の節は 13 個） | 跳び 2、2.5、3、3.5、4、4.5、5 |
| K_R0 | \([2,5]\) | 15-T1 \(=0.05\)（\(R=0\)、\(K_1>0\)） | 跳び 2、5 |

K_A3 の面積（\(Q_h=1\)）: \(A_{13a}=0.12\)、\(A_{\rm total}=0.19\)、\(\|K_{\rm total}\|_1/m=0.21\)（\([1.5,2)\) で `total` が \(-0.02\)）。

### 1.4 一次線と出力の端

| 記号 | 一次線（位置・端は eV、強度は \(I\)） |
| --- | --- |
| D1000、D1000zero、D1008、D998 | δ 1000（強度 1）、δ 1000（強度 0）、δ 1008（強度 1）、δ 998（強度 1） |
| P_A1 | δ 1000 と 1012、強度 1 と 0.5 |
| B_A1p | ビン端 999.5、1000、\(1000+2^{-10}\)、積分強度 0.6、0.4（広いビンと狭いビン） |
| B5 | ビン端 998、998.5、999、999.5、1000、1000.5、積分強度 0.1、0.2、0.4、0.2、0.1 |

| 記号 | 出力の端（eV） |
| --- | --- |
| E_A1 | 960.3、970.1、980.7、985.2、990.9、993.4、996、998.5、1000、1001（非整列。9 ビン） |
| E_A1s、E_A1s_plus8 | 988、990、994、995、996、998、1000、1001／全点に 8 を足したもの |
| E_A3 | 968、972、976、980、984、988、990、992、994、995、996、997、998、999、1000、1001（15 ビン） |
| E_A5 | 962、970、976、980、983、985、987、989、991、993、995、997、999、1000、1001 |
| E_B120 | `{"min": 880.5, "max": 1000.5, "step": 1}`（120 ビン。端は参照の格子の点に乗る。3次は低い側の外にも届く） |
| E_M | 960、980、984、988、995、999、1000、1001（M4a の E_RB1） |
| E_T1_eq、E_T1_up、E_T1_down | 997、998、999、1000、1001／第1点を `nextafter(997, +inf)`（`0x1.f280000000001p+9`）／`nextafter(997, -inf)`（`0x1.f27ffffffffffp+9`） |
| E_T2 | 950、970、980、990、999、1000、1001 |
| E_SGN_rel、E_SGN_abs | 989、990、992、993、993.125、994、996、1000、1001／989、990、992.5、993.5、996、1000、1001 |
| E_D100、E_DS2、E_DS0 | 900、1000、1001／997.5、998.5、999.5、1000.5／1000、1001、1002 |

丸めに依る入力（1.3、7.7、0.0072 など）の16進は表にある。K_T4 の3つと E_T1 の3つは二進で厳密に作る。

## 2. 既定の精度と資源の上限

| 引数 | 既定 | 出所 |
| --- | --- | --- |
| `epsrel`、`area_epsabs`、`spectrum_epsabs_fraction`、`panel_limit`、`max_cells`、`kernel_epsabs`、`kernel_epsrel`、`kernel_limit` | \(10^{-6}\)、\(10^{-10}\) Å⁻¹、\(10^{-9}\)、400、\(2\times10^6\)、\(10^{-12}\)、\(10^{-9}\)、200 | M4a のまま |
| `order_epsrel`、`order_epsabs_fraction` | \(10^{-6}\)、\(10^{-9}\) | 設計 §14 |
| `order_limit` | 64 | 設計 §14 |
| `max_convolution_points` | \(2\times10^8\)（17×交差区間の数の、2次から \(N\) 次までの累計） | 設計 §14 |
| `max_cells` | \(2\times10^6\)（1次は M4a のまま。高次は次数ごと） | 設計 §14 |
| `neglect_tolerance_fraction` | 既定値なし（`None` は判定しない） | U4 |
| `model`、`max_order` | 必須（既定値なし） | U1、U4 |

- `panel_limit` は高次の網にも次数ごとに当てる（二分で超える前に検査し、超えるなら `order_budget_not_met`）。
- **初めの網:** 必須の節だけ（§3.2）で作り、二分に任せる。
- **正式の畳み込み:** 各次数で最終の網に行う。回数を `order[n].mesh.formal_folds` に残す（数は凍結しない）。
- **和:** 合計のスペクトルの誤差は、\(J_0\) と各次数の返す誤差の絶対和に加算の丸めを足して返すだけにする。和そのものに予算は課さない（設計 §6.2 の候補を採った）。

**行ごとに変える設定**（書いていない引数は既定）:

| 行 | 設定 |
| --- | --- |
| RB1–RB8、RB1a | `epsrel` \(10^{-7}\)（設計 §14 の D2。継承の判定 (i) に余裕を持たせる） |
| A3f、AE1、T2e | `epsrel` 0.5、`area_epsabs` 1.0、`spectrum_epsabs_fraction` 1.0（点の誤差の注入を通すため） |
| N19 | `spectrum_epsabs_fraction` \(10^{-300}\)、`area_epsabs` \(10^{-6}\) |
| N20 | `area_epsabs` \(10^{-6}\) |
| N13、N15、C10、N16 | `max_convolution_points` 4500、50,000、50,000、4500 |
| N14、N30 | `max_cells` 300、`panel_limit` 9 |
| N18、N21 | 試験用の `order_panel_limit` 4、`order_min_panel_width_eV` 4.0（§4.6） |
| N31、N32、N33 | `order_epsrel` `nan`、`max_convolution_points` `True`、`order_limit` 0 |

**RB の行の設定の根拠:** 凍結の前の試走の確認で、RB の全行に `epsrel` \(10^{-7}\) を残した（2・3 次の継承の比の最大は 0.091。記録は `work/consult/m4b/trial/records/rb_inherited.json`）。これは設定を選ぶための試走の記録で、期待値ではない。試作の包絡は実装より粗いので、おそらく大きめの見積りである。

CLI の行は CONFIG に `"integration": {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200}` を明示する（M4a と同じ理由）。

## 3. 凍結する数え方

### 3.1 台の次数と状態

- 正の強度の一次線の台 \([p_{\min},p_{\max}]\)、出力の下端 \(d_{\min}\)、窓の下端 \(a\) から \(N_D=\max(0,\lfloor(p_{\max}-d_{\min})/a\rfloor)\)。入力の binary64 を厳密な有理数として計算し、等号の次数は残す。正の強度がなければ `support_order = 0`。窓の検査（M4a の `invalid_window`）は \(a\) での除算より前に行う（N27）。
- 状態: \(N\ge N_D\) なら `complete_on_display`。そうでなく許容を与え、表示のビンすべての上界が \(\le\) 許容 \(\times S_{\rm abs}\) なら `bounded_on_display`。許容がなければ `truncated`（成功）。両側の外は判定に入れない。
- 幾何から決まる \(N_D\) と、落とした量の値は別の欄にする。\(R\equiv0\) でも \(N<N_D\) なら `truncated` のままで、上界と全軸の残りが厳密に 0 になる（REF4a）。
- **`"display_support"` で \(N_D=0\):** 内部の M4a を1回実行する（`single_loss`、\(\hat r\)、段βの配送 1 回は返す）。`orders=()`、`computed=0`、表示と両側の外の合計は \(J_0\)、表示の上界は 0、全軸の残りは設計 §5 の \(N=0\) の式。整数の `max_order` は \(\ge1\) のまま（DS0）。

### 3.2 必須の節

次数 \(n\) の外側の網の必須の節は、各因子の特異点（窓の端と、解析の入口の `jumps_eV`・`kinks_eV`）の \(n\) 個の和である。有理数で作り、重複を除いてから1回だけ binary64 に丸める。物理の入口では \(ka+(n-k)b\)（\(k=0..n\)）。1次のパネルの内部の端（刻みの端）の和は必須にしない。節の列を `order[n].mesh.knots_eV` に返す。

### 3.3 打ち切りの上界

- 領域 \(D\)（下端 \(d_D\)、\(U_D=p_{\max}-d_D\)）の上界は \(S_{\rm abs}e^{\tau U_D}\Lambda_{K_1}(\tau)\Lambda_R(\tau)^N/(1-\Lambda_R(\tau))\) の、\(\tau\) の格子での最小。
- **\(\tau\) の格子:** \(\tau_i=10^{-3+i/20}\) eV⁻¹（\(i=0,1,\dots\)）のうち \(\tau_i\le40/a\) のもの。\(a=1\) eV で 93 点（最後は \(10^{1.6}\)）。\(\|R\|_1\) の上側の値が 1 未満のときだけ、\(\tau\to0\)（\(S_{\rm abs}\|K_1\|_1\|R\|_1^N/(1-\|R\|_1)\)）も候補にする。
- **\(\Lambda_f\) の上側の値:** \(\sum_{\rm パネル}\big(\int_{\rm パネル}e^{-\tau t}|\hat f|\,dt+e^{-\tau t_{\rm lo}}e_{\rm パネル}\big)\)。\(K_1\) は `total` の列、\(R\) は 13a の列、どちらも \(m\) の区間の上端を掛ける。\(e_{\rm パネル}\) はそのパネルの \(L^1\) の誤差（§4.6 の 8）。値の部分は、定数のパネルなら指数の積分を `expm1` などで安定に求めて上へ丸める。左端の重みは誤差の項にだけ掛ける。\(1-\Lambda_R>0\) を上側の値で判定できない \(\tau\) は使わない。対数で計算する。
- **ビンごとの3つの場合**（`bound`、`bound_reason`、`tau_per_eV`）: (a) そのビンに \(N\) より高い次数が届かない（\(U_D<(N+1)a\)。厳密な有理数。等号は届く側）なら `(0.0, "complete", None)`。状態が `truncated` でも、届かないビンにはこの場合を当てる。(b) 上界が求まれば `(値, None, τ)`。\(\tau\to0\) を採ったら \(\tau\) は `0.0`。正の上界が binary64 の下へあふれたら、最小の正の数へ外向きに丸める（厳密な 0 と区別する）。(c) 使える \(\tau\) がなければ `(None, "no_usable_tau", None)`、表せなければ `(None, "not_representable", None)`。
- 上界は「領域 \(D\) とそれより高い KE 側の全部に落ちる量」の累積の上界で、推定の誤差を入れた上界である（数学的な保証ではない）。
- **全軸の残り（`neglected_total_area`）:** どの分岐でも `{"value", "interval", "reason"}` を返す。設計 §5 の順で、(1) \(S=0\) か1回目の核がゼロなら 0、(2) \(R=0\) なら \(N=0\) で \(Sw\)、\(N\ge1\) で 0、(3) 繰り返しの核が非負（評価点と連続区間の最小が 0 以上）なら \(\hat r\) の区間で値・`diverges`・`undetermined`。繰り返しの核に負の値があれば、値は `None`、理由は `repeated_kernel_not_nonnegative`（行は置かない）。

### 3.4 資源の検査の3段

**課金の単位:** `max_convolution_points` は「実際に積分した、幅が正の交差区間の数 \(\times17\)」の、2次から \(N\) 次までの累計である。幅 0 の交差は数えない。次数ごとの実測を `order[n].mesh.quadrature_points`、累計を `diagnostics["resources"]["quadrature_points_total"]` に返す。

| 段 | 式 | 超えたとき |
| --- | --- | --- |
| (i) 入力だけ（段αの前、核の評価 0 回） | 次数: \(N>\) `order_limit`。求積点の下限: \(527\sum_{n=2}^{N}n>\) `max_convolution_points`。セルの下限（スペクトルの入口だけ）: ある \(n\le N\) で \(n\times(\)出力のビン\(+2)>\) `max_cells` | `MultipleLossInputError`（`order_limit_exceeded`、`work_lower_bound_exceeded`） |
| (ii) 段βの後（高次の配列の確保の前） | 見積り \(280.5\sum_{n=2}^{N}n(P_1+n-1)>\) `max_convolution_points`（\(P_1\) は1次の最終のパネルの数） | `MultipleLossConvergenceError("convolution_work_limit_exceeded")` |
| (iii) 適応の途中 | 実際の累計が上限を超える | `convolution_work_limit_reached` |

- **求積点の下限:** 次数 \(n\) の初めの網は \(n\) 枚以上で、各パネルの内部に固有の点が \(15+16=31\) 個あり、内部の点は幅が正の交差区間を 1 つ以上持つので \(31\times17\times n=527n\)。共有する端の点は数えない。設計 §14 の案 \(561\sum n^2\) も、第1稿の \(561\sum n\) も、実際の評価の数の下限ではなかった（Codex の指摘 5。1次が1枚の箱で、2次の初めの2枚は多くても \(64\times17=1088\) 点）。
- この下限は、通常の経路が全窓の畳み込みを必ず評価することを前提にする。材料がゼロ応答（強度がすべて 0）の構造的なゼロでは高次を評価しないので、下限の検査をしない。それ以外に、データを見てから評価を省く短絡は置かない。異なる必須の節が同じ binary64 に丸まるときは、網の理由（`order_mesh_unresolved`）で扱う。
- 等号（下限や見積りが上限に等しい）は通す。段βの後の見積りは設計の案のまま（下限ではなく見積り）。\(P_1\ge1\)、\(n\ge2\) で見積りは下限以上になる。
- `"display_support"` で \(N_D>\) `order_limit` のときも `order_limit_exceeded`（必要な次数を付ける）。

### 3.5 評価の数え方と段βの受け取り

- 公開の入口の「評価 0 回」は `xpsloss.xps_single_loss._call_kernel` が1回も呼ばれないこと。解析の入口では、渡した `kernel(T, side)` の呼び出しを数える。
- **1次は M4a の実行を1回だけ行う。**`_run(..., on_final_panels=None)` に受け取り口を足す。段βの `_integrate` が返したパネルを `_run` が持ち、整合の検査と警告の処理が済んだ成功の出口で `on_final_panels(result, panels)` を1回だけ呼ぶ。段α、二分の途中、分母が決まる前、例外のときは呼ばない。旧入口では既定の `None` のままで、演算の順は変えない。
- 試験は、製品の自己申告（`diagnostics["first_order"]["beta_deliveries"] == 1`）だけで済ませない。`_run` を包んで、公開の旧入口と新入口で `_call_kernel` の呼び出しの列（回数、配列、実引数）が同じこと、`on_final_panels` が受け取ったパネルの列（オブジェクト）が高次の入力（`panels_by_order[1]` と同じ値の配列）になったことを観測する。

## 4. 実装が守る約束（試験が使う面）

実装担当は参照値を見ない。この節と §5 の表だけで試験を書けるようにする。非公開の名前と欄の形は、Codex の1回目の確認の表で決めた。

### 4.1 公開の API（`xpsloss.xps_multiple_loss`。名前は `xpsloss` からも公開する）

```python
MULTIPLE_LOSS_MODEL_ID = "first_loss_13a_chain_v0"
multiple_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *,
    loss_window_eV, weight_length_A, epsilon_momentum, model, max_order,
    neglect_tolerance_fraction=None, order_limit=64, max_convolution_points=200_000_000,
    order_epsrel=1e-6, order_epsabs_fraction=1e-9,
    theta_deg=0.0, hole_charge=1.0, epsrel=1e-6, area_epsabs=1e-10, spectrum_epsabs_fraction=1e-9,
    panel_limit=400, max_cells=2_000_000, kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200,
    weight_length_kind="external_constant", extra_inverse_length_per_A=0.0, weight_length_seed_A=None,
    multiplier="weight_length", multiplier_A=None,
) -> MultipleLossSpectrumResult
multiple_loss_order_areas(material, electron_energy_eV, *, loss_window_eV, weight_length_A,
    epsilon_momentum, model, max_order, order_limit=64, max_convolution_points=200_000_000,
    order_epsrel=1e-6, order_epsabs_fraction=1e-9, theta_deg=0.0, hole_charge=1.0, epsrel=1e-6,
    area_epsabs=1e-10, panel_limit=400, kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200,
    weight_length_kind="external_constant", extra_inverse_length_per_A=0.0, weight_length_seed_A=None,
    multiplier="weight_length", multiplier_A=None,
) -> MultipleLossOrderAreasResult
```

- 公開する名前: 2つの関数、2つの結果の型、`MultipleLossInputError`、`MultipleLossConvergenceError`、`MULTIPLE_LOSS_MODEL_ID`。面積だけの入口は一次線、出力の端、表示の許容を取らない。
- `model` は `"first_loss_13a_chain_v0"` だけ。参照の模型の ID は公開の入口で `invalid_model`。`max_order` は `bool` でない整数 \(\ge1\) か `"display_support"`（面積だけの入口は整数だけ）。`neglect_tolerance_fraction` は `None` か有限の正の数（0 は拒否）。
- **高次の設定の検査（`invalid_settings`）:** `order_epsrel`、`order_epsabs_fraction` は有限の正。`order_limit`、`max_convolution_points` は `bool` でない整数 \(\ge1\)。
- M4b-1 は KE だけの署名で、`energy_axis` の引数を置かない（M4b-2 で足す。そのとき C09 は旧行を版つきで残す）。
- 既存の関数、型、例外、CSV、`run.json`、既定の数値は変えない。1次は `single_loss_spectrum`（面積だけの入口は `loss_window_area`）と同じ引数の M4a の実行で、検査の順と例外もそのまま。M4b の入力の検査（`model`、`max_order`、許容、高次の設定、§3.4 の (i)）は M4a の実行の前に行うが、窓の検査は \(N_D\) の計算より前に M4a の規則で行う。

### 4.2 結果の欄

**`MultipleLossSpectrumResult`**（\(n_b\) は出力のビンの数。配列は `float64`、形 \((n_b,)\)、出力の端の昇順。表のキー `bin[i]` は配列の添字 \(i\)）:

| 欄 | 内容 |
| --- | --- |
| `single_loss` | 同じ引数の `SingleLossSpectrumResult` そのもの（全欄がビット単位で同じ。旧版の `comparison` もそのまま） |
| `model_id`、`depends_on_model_id` | `MULTIPLE_LOSS_MODEL_ID`、`single_loss.model_id` |
| `kernel_energy_policy`、`repeated_kernel`、`hole_treatment` | `"fixed_E0"`、`"term_13a"`、`"first_event_only"` |
| `orders` | タプル。要素の `order` は 1..\(N\)（\(N_D=0\) の `"display_support"` では空）。表のキー `order[n]` は `order == n` の要素 |
| `zero_loss`、`zero_loss_error` | `single_loss` の同名の配列の写し |
| `spectrum`、`spectrum_error`、`outside[side]["spectrum"]`、`outside[side]["spectrum_error"]` | \(J_0+\sum_{n\le N}J_{n,\rm total}\) と、返す誤差の絶対和＋加算の丸め。`side` は `"low"`、`"high"` |
| `order_area_ratio` | `{"value", "error", "interval": (下端, 上端)}`。`value` は \(\hat mA^\beta_{13a}\) を1回丸めた値。区間は、\([m_-,m_+]\) と \([A_{13a}-e,A_{13a}+e]\) の四隅を厳密な有理数で掛け、最小を下へ、最大を上へ、最後に1回だけ丸める。`error` は中心から両端までの大きい方（上へ丸める） |
| `convention_ratio` | `"exactly_one"`、`"below_one"`、`"not_determined_by_convention"` |
| `truncation` | 下 |
| `diagnostics`、`warnings` | §4.4 |

**`orders` の要素**（`<c>` は `total`、`no_hole`、`hole`。`part` は `inherited`、`propagated`、`local`、`multiplier_interval`。面積の内訳だけ、5つ目の `convention` を持つ）:

| 欄 | 内容 |
| --- | --- |
| `order` | 次数 |
| `<c>`、`<c>_error` | ビンごとの積分強度（\(I\)）と返す誤差。負を切らない |
| `<c>_error_parts[part]` | 形 \((n_b,)\) の配列。添字は `<c>_error_parts[part][i]`。前の3つは固定の乗数 \(\hat m\) での内訳で、和が固定の乗数の誤差 \(e_f\)。`multiplier_interval` は \(\rho_n(|\hat y|+e_f)\) |
| `outside[side]` | 辞書。`<c>`、`<c>_error`、`<c>_error_parts`（`outside[side]["total_error_parts"][part]` はスカラー） |
| `area`、`area_error` | 列の名前から float への辞書。\(w_{n,j}\)（無次元）と返す誤差 |
| `area_error_parts[part][column]` | 辞書の辞書。部分が外、列が内。5つの部分の和（上へ丸める）が `area_error`。`convention` は、規約から来る旧版の誤差を1回だけ載せる欄で、S2 の次数 1 の `total` のほかは 0 |
| `area_zero_kind[column]` | `None`、`"exact_zero"`、`"rounded_zero"` |
| `cutoff_sensitivity`、`cutoff_sensitivity_error` | 列の辞書。\(h_{n,j}\)（無次元。固定の \(m\)、\(L_w\) の偏微分） |
| `mesh` | `knots_eV`（必須の節）、`initial_panel_edges_eV`、`panel_edges_eV`（最終の実際の網）、`panels`、`stages`（二分の段数）、`history`（段ごとの枚数と、正式か予備かの判定の結果）、`builds`（その次数の網を作り始めた回数。1 のはず）、`local_area_terms`（最終の網のパネルごとの局所の面積の誤差。11 列。パネルの単位）、`formal_folds`、`quadrature_points`（その次数の実測）、`cells`。面積だけの入口では `cells` と `formal_folds` は 0 |

- **固定の乗数と返す誤差:** パネルは \(m\) を掛けない値で持ち、有理数の \(\hat m^n\) を最後に1回掛ける。2次以降は固定の乗数の誤差だけを再帰に渡し、各次数の最終の出力に1回だけ \(\rho_n=\max(1-(m_-/\hat m)^n,(m_+/\hat m)^n-1)\)（厳密な比の累乗から上へ丸める）を当てて、返す誤差を \(e_f+\rho_n(|\hat y|+e_f)\)（上へ丸める）とする。乗数に区間がなければ \(\rho_n=0\)。
- **次数 1 は全次数の式から外す（旧欄の転記）:**
  - ビンと両側の外の値と返す誤差は、`single_loss` の `single_loss_<c>`、`single_loss_<c>_error` をビット単位で写す。\(N=1\) の `spectrum`、`spectrum_error` と両側の外も、旧版の欄（旧版の丸めの順）をそのまま写す。
  - 内訳: `inherited` は区間の補正の前の `single_loss.diagnostics["error_budget"]["outputs"][key]["estimate"]`、`propagated` と `local` は 0、`multiplier_interval` は「返す誤差 − `inherited`」を上へ丸めた値（区間がなければ 0）。
  - 面積の `total`: `area["total"]` は `single_loss.single_loss_weight`、`area_error["total"]` は `single_loss.single_loss_weight_error`（ビット単位）。内訳は、乗数に区間がなければ `inherited` が全部。区間があって S2 でなければ `inherited` が固定の乗数の誤差、`multiplier_interval` が残り。**S2 では**旧版の規則で返す誤差が \(|\hat w_1-1|\) の上向きの丸めになり（A3f では 0）、固定の乗数の面積の誤差 \(\hat me_A\)（A3f で約 0.0368）より小さいことがある。差を `multiplier_interval` に入れると負になるので、S2 の次数 1 の `total` だけは `convention` に旧版の誤差を載せ、ほかの4つの部分を 0 とする（\(w_1=1\) は規約による制約。`diagnostics["first_order"]["area_constraint"]` に `"inverse_window_area"` と書く）。
  - **固定の乗数の再帰の種は別に供給する:** `diagnostics["first_order"]["fixed_multiplier_area_error"][column]`（11 列）は、段βのパネルの面積の誤差から \(\hat me_{A_j}+\epsilon|\hat mA_j|\) を上へ丸めた値で、返す誤差の規約（S2 の \(|\hat w_1-1|\)、区間の補正）を含まない。2次以降の Young の式と \(\hat r\) の固定の誤差はこれを使う。試験は、この欄（`total` を含む通常の式）と、上の S2 の例外を別々に照合する。
  - 面積のほかの 10 列（旧版に無次元の欄がない）: `area[j]` は \(\hat mA_j\) を厳密な有理数で掛けて1回丸めた値。固定の乗数の誤差は \(\hat me_{A_j}+\epsilon|{\rm area}_j|\) を上へ丸めた値（`inherited`）。返す誤差は旧版の \(\rho\)（`single_loss.diagnostics["multiplier"]["rho"]`。区間がなければ 0）で \(e_f+\rho(|{\rm area}_j|+e_f)\) を上へ丸めた値。
  - 感度 \(h_{1,j}\) と誤差は `single_loss.cutoff_sensitivity`、`cutoff_sensitivity_error` の写し。
- **次数 \(n\ge2\) の面積と感度:** \(w_{n,j}\) はその次数のパネルの面積（\(\hat m^n\) を掛けたもの）。恒等値 \(w_{1,j}\hat r^{n-1}\) で上書きしない。\(h_{n,j}=h_{1,j}\hat r^{n-1}+(n-1)\,w_{1,j}\,h_r\,\hat r^{n-2}\)（\(h_r=h_{1,13a}\)。\(\hat r=0\) で負のべきを評価しない形）。
- **高次の感度の誤差（`cutoff_sensitivity_error`）:** 入力は区間の補正の**前**の、固定の乗数 \(\hat m\) での4つの量とその誤差である。\(h_{1,j}\) と \(h_r\) は窓の下端の節の値と点の誤差から（\(\hat m\,a\,e_{K_j}(a)+2\epsilon|h|\) を上へ丸める）、\(w_{1,j}\) は上の再帰の種、\(\hat r\) は \(\hat mA_{13a}\) と、同じ再帰の種 `fixed_multiplier_area_error["term_13a"]`（\(\epsilon|\hat mA_{13a}|\) の項を含む）。感度と再帰は同じ丸めの規約を使う。この4つの区間で上の式の和と積を外向きに評価し、中心からの最大の距離を固定の乗数の誤差 \(e_f\) とする。最後に乗数の区間を1回だけ足して \(e_f+\rho_n(|h_{n,j}|+e_f)\)（上へ丸める）を返す。旧版が区間を補正した後の1次の感度の誤差（S2 など）を、この計算に入れ直さない。

**`truncation`:** `requested`（整数か `"display_support"`）、`computed`、`support_order`、`state`、`neglect_tolerance_fraction`、`neglected_at_or_above_region_bound`（長さ \(n_b\) のリスト。`None` を許す）、`bound_reason`（同じ長さ）、`tau_per_eV`（同じ長さ）。3つの欄の組は §3.3 の (a)–(c)。`required_order_estimate`（許容を与えたとき。満たす最小の次数の見積り。求まらなければ `None`）、`neglected_total_area`（`{"value", "interval", "reason"}`。§3.3）。

**`MultipleLossOrderAreasResult`:** `single_loss`（`LossWindowAreaResult`）、`model_id`、`depends_on_model_id`、`kernel_energy_policy`、`repeated_kernel`、`hole_treatment`、`orders`（`order`、`area*`、`area_zero_kind`、`cutoff_sensitivity*`、`mesh`）、`order_area_ratio`、`convention_ratio`、`diagnostics`、`warnings`。`truncation` は持たない。

- 一次線の強度がすべて 0 なら、ビンと外は厳密に 0（誤差も 0）、面積は計算する（CT1）。
- 材料がゼロ応答なら M4a の構造的なゼロの経路で、高次の値・誤差・面積は厳密に 0、\(J=J_0\)（M4A-Z）。S2 との組は M4a-2 の `multiplier_undefined`、窓の 13a の重みで \(c=0\) なら `weight_length_undefined` が先（N28、N29）。
- **`area_zero_kind`:** `"exact_zero"` は、構造を示す証拠（強度がゼロの材料、恒等的にゼロと分かる繰り返しの因子、一次線に依らない構造的なゼロ）がある場合だけ。任意の呼び出し可能な核の数点の値から、恒等的なゼロや反対称を認定しない。`"rounded_zero"` は、中心の値が binary64 の 0 で証拠がない場合。ほかは `None`（A3z の `hole` は入力を二進の有理数として読むと約 \(3.5\times10^{-18}\) で、`None` か `"rounded_zero"`。床の内であることだけを見る）。
- 解析の入口では、\(E_0\) に関わる欄（`single_loss.electron_energy_eV`、`diagnostics["cumulative_loss_over_E0"]`）は `None`。

### 4.3 例外、CLI、記録

| 型 | 状態 | 理由（`reason`） |
| --- | --- | --- |
| `MultipleLossInputError(ValueError)` | `error` | `invalid_model`、`invalid_max_order`、`invalid_neglect_tolerance`、`invalid_settings`、`order_limit_exceeded`、`work_lower_bound_exceeded`、`invalid_schema`（CLI） |
| `MultipleLossConvergenceError(NumericalConvergenceError)` | `not_converged` | `order_budget_not_met`、`inherited_error_exceeds_order_budget`、`order_mesh_unresolved`、`convolution_work_limit_exceeded`、`convolution_work_limit_reached`、`truncation_bound_exceeds_tolerance`、`truncation_bound_unavailable`、`nonfinite_output` |
| M4a・M4a-2・上流の例外 | そのまま | 型、`reason`、`phase`、`stage`、`inner_reason`、公開の属性を保つ（N27–N30） |

**属性**（該当しないものは `None`）:

| 属性 | 持つ理由 | 内容 |
| --- | --- | --- |
| `reason` | すべて | 上の文字列 |
| `argument` | `invalid_settings` | 引数の名前 |
| `requested_order`、`order_limit` | `order_limit_exceeded` | 要求した（または必要な）次数と上限 |
| `resource`、`lower_bound`、`limit`、`order` | `work_lower_bound_exceeded` | `"max_convolution_points"` か `"max_cells"`、下限、上限、セルの場合は最初に超える次数 |
| `order`、`stage`、`inner` | 収束の例外のすべて | 次数（次数に属さなければ `None`）、段、内側の原因。`stage` は `"after_first_order"`、`"order_refinement"`、`"order_fold"`、`"truncation"`。`inner` は `None` か、`{"type", "reason", "stage"}` を持つ JSON にできる辞書 |
| `estimate`、`measured`、`limit` | 資源の2つ | 段βの後の見積り（`convolution_work_limit_exceeded`）、止まった時点の実測（`convolution_work_limit_reached`）、上限 |
| `exceeded`、`key`、`side`、`column`、`suggestion` | `inherited_error_exceeds_order_budget` | `exceeded` は超えた全キーのリスト。各要素は `key`、`inherited`、`budget_relative`、`budget_absolute`、`side`、`suggestion`。代表（`key` ほか）は比が最大のもので、同率なら面積の列の順、領域の順、成分の順（`total`、`no_hole`、`hole`）で先のもの |
| `required_order_estimate` | `truncation_bound_exceeds_tolerance` | 見積り |
| `panel_eV`、`panel_limit` | `order_mesh_unresolved`、`order_budget_not_met` | 二分できなかったパネルの両端、当たった枚数の上限 |
| `diagnostics` | 収束の例外 | 止まる直前のスナップショット |

- **(i) の側:** 予算は \(\max(B_{\rm abs},B_{\rm rel})\)。`side` は大きい方の側（同じなら `relative`）。`suggestion` は、相対の側が `"epsrel"`、絶対の側がビンなら `"spectrum_epsabs_fraction"`、面積なら `"area_epsabs"`。文言は、1次の網と核の絶対の精度を確かめる案内にし、予算を縮めるだけで直るとは書かない。
- **止まる理由の使い分け:** (iii) を満たさないまま、次の二分で高次の枚数の上限（`panel_limit`、試験では `order_panel_limit`）を超えるなら `order_budget_not_met`。二分が要るのに子の幅が表せない（節点の分離が成り立たない、または試験用の `order_min_panel_width_eV` より狭い）なら `order_mesh_unresolved`。
- **`nonfinite_output`:** 最終の出力（値か誤差）が有限でないときに出す。大きな \(\hat m^n\) の途中の丸めだけでは出さない（係数は有理数で持つ）。
- 打ち切りの判定は、要求した次数までの計算と判定が済んだ後に行う（`stage` は `"truncation"`、核の評価は 1 回以上）。

**CLI:**

```powershell
xpsloss xps-multiple-loss CONFIG --spectrum SPECTRUM.json --loss-window-eV A B --weight-length-A L --epsilon-momentum M --model first_loss_13a_chain_v0 --max-order N --out DIR
```

- `--model`、`--max-order`（整数か `display_support`）は必須。任意: `--neglect-tolerance-fraction`、`--order-limit`、`--max-convolution-points`、`--order-epsrel`、`--order-epsabs-fraction`。ほかの引数と `spectrum.json` は `xps-single-loss` と同じ（`axis` は `"kinetic_energy_eV"` だけ。ほかは `invalid_schema`）。
- 出力: `multiple_loss_spectrum.csv`（見出しは `region,edge_low_eV,edge_high_eV,zero_loss,zero_loss_error`、続いて次数ごとに `order<n>_total,order<n>_total_error,order<n>_no_hole,order<n>_no_hole_error,order<n>_hole,order<n>_hole_error`、最後に `spectrum,spectrum_error`。行は `low`、`0`…（ビンの昇順）、`high`。外の無限の側の端は空欄）、`order_areas.csv`（見出し `order,column,area,area_error,cutoff_sensitivity,cutoff_sensitivity_error`。行は次数、その中で 11 列の順）、`run.json`。浮動小数は `%.17g`（読み直して、符号つきのゼロまで照合する）、ASCII。全 CSV の後に `run.json` を書く。
- 成功は終了 0。失敗は終了 2 で、出力先を作らず、stderr の JSON に `status`、`error_type`、`reason`、`message` と上の属性のすべてを残す。既存の出力先は `FileExistsError` で拒否し、既存のファイルのバイトを変えない（C07）。構文エラー（必須の引数の欠落）は M4a と同じく `ValueError` で `reason` なし。
- **`run.json`:** 既存の共通の欄に加えて、`model_id`、`depends_on_model_id`、`dependent_model_ids`（M4a の3つのキーに `single_loss` を足した辞書）、`repeated_kernel`、`hole_treatment`、`kernel_energy_policy`、`orders_definition`（次数ごとに、次数、厳密な係数 \(\hat m^n\)（整数の比の文字列）、乗数、区間）、`truncation`、`order_area_ratio`、`convention_ratio`、`effective_settings` と `effective_settings_sources`（M4a の 8 個に `order_epsrel`、`order_epsabs_fraction`、`order_limit`、`max_convolution_points`、`neglect_tolerance_fraction` を足した 13 個。出所は `cli`、`config`、`default`）、`spectrum_definition`（M4a-2 の規約の欄に `max_order` を足す）、`results`（§4.2 と同じ名前）、`comparison`（群 `orders`、`areas`、`spectrum` がすべて `not_matched`。入れ子の `results.single_loss` の `comparison` は旧版の結果のまま）、`physical_validity: "not_established"`、`runtime`、`warnings`、入力（設定、スペクトル）の SHA-256。BE の変換の欄は M4b-2。
- **再生（C01）:** `run.json` の `config`、2つの入力の SHA-256、`spectrum_definition`、出力の端、`max_order`、`effective_settings` から API を再生し、値、誤差、状態、警告、模型の ID、単位、版を照合する。
- **`capabilities`:** `supported` は既存の要素と順のまま、末尾に `"xps_multiple_loss"` を足す。基準のコードにあるブロック（サブコマンドごとの `unsupported_physics` を含む。最上位の `bulk_unsupported_physics` はバルク用）は基準と等しいまま（`version` は §5.2 の規則）。新しいブロック `multiple_scattering` は `model_id`、`emission: "normal_only"`、`kernel_energy_policy: "fixed_E0"`、`status: "approximate"`、`unsupported_physics`（`oblique_emission`、`trajectory_dependent_weight_length`、`self_consistent_length`、`intrinsic_loss_separation`、`losses_outside_window`、`elastic_scattering`、`instrument_broadening`、`background_subtraction`、`energy_dependent_kernel`、`binding_energy_axis`）を持つ。参照の模型の ID はどこにも出さない。

### 4.4 `diagnostics` と警告

`diagnostics` は、`first_order`（`beta_deliveries`、`final_panels`、`area_constraint`、M4a の `kernel_calls` の写し）、`resources`（`work_lower_bound_points`、`work_estimate_after_first_order`、`quadrature_points_total`、`cells_max`）、`negative_values`（`total` と 13a の核の、評価点の最小と連続区間の最小。連続区間の最小は、パネルの端と、補間多項式の内部の極値を含めて求める）、`cumulative_loss_over_E0`、`effective_settings`、`warnings` を持つ。

外側の `warnings` は文言のタプル。`diagnostics["warnings"]` は `{"source", "code", "message"}` の、順序のあるリスト。試験は識別子で判定し、常に出るものと条件つきのものの集合を照合する。

| 識別子（出所 `xps_multiple_loss`） | 条件 |
| --- | --- |
| `model_is_approximation`、`excess_enters_once`、`depth_correlation_dropped`、`no_elastic_scattering`、`order_areas_not_probabilities` | 常に |
| `truncated` | `truncation.state` が `truncated` |
| `order_area_ratio_not_below_one` | \(\hat r\) の区間の上端が 1 以上 |
| `multiplier_differs_from_weight_length` | \(m\ne L_w\)（binary64） |
| `cumulative_loss_over_E0_gt_0p1` | \(Nb/E_0>0.1\)（公開の入口だけ） |
| `negative_total_kernel`、`negative_13a_kernel` | 評価点か連続区間の最小が負 |
| `dispersive_k_parallel_13a` | `k_parallel` で、強度が正の振動子に分散 \(\ne0\) のものがある |
| `E0_above_5keV` | \(E_0>5000\) eV |

M4a・M4a-2 から継承した警告（`window_only`、`kernel_fixed_at_E0` など）は、旧い出所（`xps_single_loss`、`xps_depth_averaged`）のまま持つ。

### 4.5 参照の模型（非公開）

`xpsloss._loss_reference_models.integrate_orders_reference(kernel, loss_window_eV, *, model, ...)`。引数は §4.6 の 1 と同じで、`model` は `full_kernel_geometric_ref`（\(C_n=K_1^{*n}\)）、`full_kernel_factorial_ref`（\(K_1^{*n}/n!\)）、`excess_poisson_ref`（\(\sum_nC_nz^n=e^{z\bar g}/(1-zR)\)）。

- **成分の規約（3つの模型に共通）:** 模型を `total` の核と `no_hole` の核にそれぞれ適用し、`hole` は差 `total − no_hole` と定義する（\(K_1\) に線形でないため）。**この非公開の入口は \(Q_h=1\) だけを受け付け、ほかの正孔電荷は `ValueError` で拒否する**（\(Q_h=0\) では単位電荷の `hole` が 0/0 になる）。B1 の `hole` は \(K_1\) に線形で、\(Q_h=0\) でも定義される（RB8、A3b）。
- 結果は `reference_only` が `True` で、`orders` の `<c>`、`<c>_error`、`outside`、`area`（3成分だけ）、`truncation` を持つ。8項の分解、B1 の内訳の3分け、B1 専用の面積の比（`order_area_ratio`、`convention_ratio`）は持たなくてよい（試験は求めない）。
- 打ち切りの上界は模型ごと（設計 §4.5。幾何は次数の因子 \(\Lambda_{K_1}^n\)、階乗は \(\Lambda_{K_1}^n/n!\)）。API、CLI、`capabilities`、`xpsloss` の公開の名前に出さない。

### 4.6 試験用の口

1. **`xpsloss._loss_orders.integrate_orders(kernel, loss_window_eV, *, features, weight_length_A, hole_charge=1.0, primary=None, output_edges_eV=None, model, max_order, neglect_tolerance_fraction=None, order_limit=64, max_convolution_points=200_000_000, epsrel, area_epsabs, spectrum_epsabs_fraction, panel_limit, max_cells, order_epsrel=1e-6, order_epsabs_fraction=1e-9, multiplier="weight_length", multiplier_A=None, initial_edges=None, inject=None, order_inject=None, order_panel_limit=None, order_min_panel_width_eV=None) -> (result, panels_by_order)`。**解析の入口。既存の `integrate_window` と同じ `kernel`、`features`、`initial_edges`、`inject` を受け、公開の経路と同じ高次の駆動を呼ぶ。繰り返しの核は同じ `kernel` の `term_13a` の列。`primary` が `None` なら面積だけ。`panels_by_order[n]`（\(n=1..N\)）はパネルのリストで、左端の順、\(m\) を掛けない値、各 `values.shape == (11, 17)`、列の順は 11 列の順。
2. **`inject(role, T, values, errors)`**（M4a と同じ）で入れた点の誤差は1次のパネルの誤差になり、高次では継承の分 (i) になる。**現行の M4a は、検証点の誤差を検証の項（幅 \(\times\) 最大）に、補間点の誤差を Clenshaw–Curtis の重みの項に別々に足す**ので、両方の役割に同じ定数 \(e\) を入れると面積の誤差は \(2(b-a)e\) になる。面積の誤差を \((b-a)e\) にしたい行（A3f、AE1、T2e）は、M4a の N27 と同じく補間点だけに入れ、検証点には元の値と誤差を返す。
3. **`order_inject(order, role, panel_eV, nodes_eV, values, local_errors) -> (values, local_errors)`。**次数 \(n\ge2\) のパネルごとに呼ぶ。`role` は `"interpolation"`（`float64[11, 17]`）か `"verification"`（`float64[11, 16]`）、`panel_eV` はパネルの（左端, 右端）、単位は \(m\) を掛けないパネルの単位。返した `local_errors` は局所の分 (iii) に入る。検証点の値を変えると局所の残差になる。
4. **`order_panel_limit`**（既定は `panel_limit`。高次の各次数の最終の枚数の上限）と **`order_min_panel_width_eV`**（既定は `None`。子のパネルの幅の下限）。試験専用のキーワード。
5. **`xpsloss._loss_orders.convolve_panels(first_panels, loss_window_eV, *, model, max_order, features, order_epsrel=1e-6, order_epsabs_fraction=1e-9, order_panel_limit=400) -> panels_by_order`。**層C。与えた1次のパネル（`Panel.from_lobatto`）から高次のパネルを作るだけで、核も一次線も使わない。乗数は 1。A4a・A4c は `features` に窓の端の跳びだけを渡す。
6. **`xpsloss._loss_orders._fold_order(panels, geometry, *, multiplier_power, error_parts) -> folded`。**高次用の領域の積分。既存の `Geometry`（有理数の領域の境と質量）と 9 点のセルを共有し、衛星の3成分を直接に積分して、固定の乗数での3つの内訳だけを返す。\(\hat m^n\) は有理数で渡す。`error_parts` は `panels` と同じ長さのリストで、各要素は `{"inherited", "propagated", "local"}` から `float64[11, 17]`（節の誤差。\(m\) を掛けないパネルの単位）への辞書。試験はこの関数を包んで、渡った節の内訳を直接に読む（A3f の行内の対照）。`zero_loss`、`spectrum`、感度のキーは返さない。検証の残差は局所の分だけ、節の誤差は出所ごとに1回だけ加える。面積の Young の項をビンに足さない。
7. **正式の畳み込みと判定の順:** 二分の候補は安い見積りで選び、最終の候補の網で正式に畳んで (i) と (iii) を判定する。(iii) を超えたときだけその次数を細分して畳み直す。(i) を超えたら止める。(ii) は実額を返す。`formal_folds` は正式に畳んだ回数で、成功したスペクトルの各高次で 1 以上、面積だけの入口で 0。前の次数の網には戻らない。
8. **包絡と \(L^1\) の誤差の種（推定であって保証ではない）:**
   - 1次のパネルの点の誤差の包絡（列ごとの定数）: \(\epsilon_{1}=\max(\max_i{\tt errors}_i,\ {\tt epsilon})+{\tt rho}+{\tt coeff\_diff}\)（`Panel` の欄。節と検証点の不確かさの最大、検証の残差、補間の係数の差）。
   - 1次のパネルの面積の誤差 \(e_{\rm パネル}\): 現行の `_areas` がそのパネルに置く5項（`interpolation`、幅\(\times\)`verification`、`errors` の Clenshaw–Curtis の重みつきの和、0、丸め）の和。乗数を掛ける前の Å⁻¹ の値を、Young の式と打ち切りの上界に使う。点の包絡に幅を掛けて作り直さない（両方の役割に誤差を入れた行で別の量になる）。**全パネルの和が `window_area_error` に等しいのは、実数の式としてである。**旧版は5項のそれぞれを先にパネルの方向に足し、最後に `math.fsum` する。パネルごとに5項を足してから全パネルを足すと、列によって 1 ulp 違う（Codex の実測）。旧版の欄とのビット単位の照合には旧版の集計（旧版の欄の写し）を使い、パネルごとの種の和を使わない。この5項は非負の誤差の推定で、任意の符号の変化に対する厳密な \(L^1\) の上界ではない。
   - 次数 \(n\ge2\): 幅が正の各交差区間で、前の次数のパネルの定数の包絡を使い、新しい網の 17 の補間点と 16 の検証点で設計 §6.2 の積分を成分ごとに求め、その最大を新しいパネルの包絡にする。\(|\hat C_{n-1}|\epsilon_R\) と混合の項は (i) に1回、\(\epsilon_{n-1}|\hat R|\) は前の次数の (i)/(ii)/(iii) を今の次数の (i)/(ii)/(ii) へ。セルの重み（補間の基底の絶対値）はビンへ移すときだけに使い、次の次数へ運ばない。
   - **面積の誤差の推定量は Young の右辺そのもの**（外向きに丸める）: \(e_n=\|\hat R\|_1e_{n-1}+\|\hat C_{n-1}\|_1e_R+e_{n-1}e_R+e^{\rm local}_n\) を、点の包絡と同じ規則で (i)/(ii)/(iii) に配る。\(\|\cdot\|_1\) は、区分定数と区分1次式なら符号の零点で分けて端と原始関数で求める。一般の多項式には上側の評価を許す。
9. 試験が `panels_by_order` から作る量（A2、A3z）: 各パネルの 17 節の値の補間多項式を、試験の側の 24 点のガウス則で積分して \(\int e^{-\tau s}\hat C_n\,ds\)、\(\int\hat C_n\,ds\)、\(\int s\hat C_n\,ds\) を求める（製品の畳み込みを使わない）。対象の列は A2 が `term_13a`、A3z が `hole`。\(m\) を掛けないパネルの値を、\(m\) を掛けない閉形式と比べる（物理の単位の出力と混ぜない）。試験の側の求積の誤差（点の数を変えた差）を \(u_r\) に足す。

## 5. 期待する状態（109 行）

**表:** `benchmarks/numerics/m4b1/expected_states_m4b1.json`（スキーマ `xpsloss-m4b1-expected-states-v1-draft4`、ASCII、LF。SHA-256 は凍結の時に統括担当が付録に書く）。生成は `python benchmarks/numerics/m4b1/expected_states_m4b1.py`（標準ライブラリだけ。`xpsloss` を読み込まない）。`--check` で再生成してバイト単位で比べる。生成は、行の総数 109、CI 90、全件だけ 19 を、登録した定数と照合する。

### 5.1 定義と数

- **行**は「新しい入口」（入口、入力、設定の組が他の行と違い、自分の期待する終了状態を持つ実行）を1つ持つ。**行内の検査**（`in_row_checks`）と**再利用**（`reuse`。`RB1.result` のように成果物まで書く）は行に数えない。CI の行は全件だけの行を再利用しない。
- **数えないもの:** skip、比較 0 件、期待と違う例外は、合格にも検出にも数えない。
- **キーの台帳（`check_plan`）:** 行ごとに、予定するキーの群を `{"keys", "kind", "oracle", "count"}` で持つ。予定したキーは、結果が `no_reference`、`unresolved`、`not_run` でも台帳に残す。

| 群 | 行 | 入口 | 行数 | CI | 分母 |
| --- | --- | --- | ---: | ---: | --- |
| RM-A1 | A1、A1s、A1p | 解析 | 3 | 3 | 数値（丸めの水準） |
| RM-A2 | A2 | 解析 | 1 | 1 | 数値 |
| RM-A3 | A3a–A3f、A3z、AE1、OI1、OI2 | 解析（A3f、AE1 は面積だけ） | 10 | 10 | 数値と誤差の内訳 |
| RM-A4 | A4a、A4c（層C）、A4b | 層C、解析 | 3 | 3 | 数値 |
| RM-A5 | A5 | 解析 | 1 | 0 | 数値 |
| RM-B | RB1–RB8、RB1a | 公開 | 9 | 1（RB1） | 数値 |
| I | I-shift、I-areas | 解析 | 2 | 2 | ビット一致、恒等式 |
| INV | INV1、INV2 | 解析 | 2 | 2 | 数値 |
| REF | REF1–REF3、REF4a–c | 参照の模型（REF4a は解析） | 6 | 6 | 数値 |
| M4A | M4A-c1–c7、-o3、-o3s、-Al、-Z、-area-c1、-area-c4 | 公開 | 13 | 7 | ビット一致 |
| T | T1a–c、T2a–c、T2e、T4u、T4d、T4v、T5、T5n、DS0、DS1 | 解析（DS0 は公開） | 14 | 14 | 状態、数値 |
| N | N01–N16、N18–N21、N26–N33 | 公開、解析 | 28 | 28 | 状態 |
| C | C01–C11 | CLI | 11 | 11 | 再生と状態 |
| CT（対照） | CT1、CT3–CT7 | 解析、公開 | 6 | 2 | 数値、状態 |
| 計 | | | **109** | **90** | |

- 期待する終了状態の内訳: 成功 68、CLI の完了 3、`error` 25、`not_converged` 13。
- **予定のキー（表の `counts`）。種類ごとに別の分母で数える:**
  - 外部の参照を読むキー 9,668: 数値の比較 `direct_reference` 9,608、上界 `bound_reference` 38、誤差の内訳 `error_part` 6、区間 `interval_contains` 16。群ごとには RM-B 6,096、RM-A 2,597、REF 534、T 253、I 88、CT 66、INV 34。違うキーが同じ参照の量を読むことがあるので、独立な条件の数ではない。
  - 報告した値からの作り直し `recomputed_from_report` 40（A3f の \(\rho_n\) の 22、1次の写像と区間の両端の 13、3次の3つの内訳の 3、AE1 の局所の分の 1、DS0 の 1）、状態と内訳の不等式 `inequality` 8（OI1 の 3、OI2 の 4、AE1 の予算の 1）。A3f の3次の作り直しは、報告した2次の内訳と、参照の \(\|R\|_1\)、\(\|\hat C_2\|_1\) を使う。
  - 恒等式 5,408、厳密な 0・構造 341、ビット一致 37、再生 3、診断 79、終了状態 109。
- **参照の記録の書式:** 参照のファイルは、`(行, キー, 種類, オラクル)` ごとに、中心の値と \(u_r\)（十進の文字列）、区間の端、比較幅、「値なし」の理由を持つ（表の `reference_record_example` が値のない見本）。参照の担当は、凍結した表のキーの群を具体のキーに展開し、展開の数が表の `count` と合うことを確かめる。精度の比較（§7.2）と推定の検査（§7.3）は同じ記録を一度読み、別々の判定の欄に残す。
- 恒等式のうち、設計 §10 の I 群の「面積の積」「全ビンと両側の外の和」「成分の関係」「台の外の厳密な 0」「ラプラス変換」「モーメント」は、A1–A5 と RM-B の行の行内の検査とキーで持つ（行を足さない）。結合則は参照の生成器の自己検査として記録する（§6.1）。A1p は有限の幅の一次線を参照と直接に比べる行で、「δ への極限」は主張しない。
- **キーの群の展開（生成のときに行い、表には並べない）:** 生成は、キーの群を具体のキーに展開して検査する。(1) 空の展開と、行の中の `(キー, 種類, オラクル)` の重複を拒む。(2) 既定値を展開し、説明だけの欄を除いた、型つきの入口の署名で重複を見る（§4.1 の既定値をすべて展開するので、省略と、明示の `multiplier="weight_length"`・`multiplier_A=None` は同じ署名になる。本物の c2（`external`、20 Å）や、真偽と整数の違いは残る。CLI の行は生の引数のまま）。(3) 再利用する成果物（`.result`、`.baseline`）を、再利用元の行が実際に供給することを見る（表の `provides`）。(4) 各変異を、行の台帳にある具体の `(行, キー, 種類, オラクル)` か、行の終了状態に結ぶ。台帳にないキー（次数や領域が行の外）と重複を拒む。(5) CI の行を回した各段階（90／19、89／20、88／21、87／22）で、件数と再利用の向きを見る。(6) 台帳の全部のキーと変異の全部の検出先を、その行の `max_order` と出力の端から決まる範囲（次数、ビン、領域）に照合する。行内の対照のキーは、表の `controls` に書いた対照の入力で照合する。(7) 自己検査として、`TableError` になるはずの反例を生成のたびに実行する: c2 に既定の乗数を明示したもの（c1 と重複）、A1s の入力だけを \(N=2\) にして3次のキーを残したもの、台帳と X09 の両方に 99 次のキーを置いたもの、出力の端を越えるビンの添字。外部の参照を読むキーは、すべて意味のある名前に展開する（A3z のモーメント 4、REF3 の 17 領域×2成分、REF4b・REF4c の各 7 ビンを、通し番号から直した。件数は同じ）。文章だけの群（診断の一部など）は、添字で数えるだけである。

### 5.2 行の内容

**解析の核（\(L_w=m=2\) Å、\(Q_h=1\)。変えた点だけ書く）**

| ID | 核、次数、一次線、端 | 期待（成功の行は `truncation.state` つき） |
| --- | --- | --- |
| A1 | K_A1、\(N=6\)、P_A1、E_A1 | 成功、`truncated`、`support_order` 39。区間ごとに多項式なので丸めの水準で比べる（§7.2）。両側の外が非ゼロ（1次は高い側の外、6次は低い側の外に届く）。`knots_eV` が有理数の節を1回丸めた列。3次の \({\rm RN}(2a+b)={\rm RN}(10.3)\) が、初めの網と最終の網の実際のパネルの端にビット単位である |
| A1s | K_A1s、\(N=3\)、D1000、E_A1s | 成功、`truncated`、`support_order` 6。`bin[2]`（損失 5–6 eV）は2次だけが入る。各次数で、台 \([na,nb]\) と正の幅で重ならない領域（1・2・3 次で 7、6、6 個）は、3成分の値も誤差も厳密に 0（57 組）。ゼロ損失のビン `bin[6]` では \(J_0=1\)、2次は厳密に 0、`spectrum` は 1 |
| A1p | K_A1s、\(N=3\)、B_A1p、E_A1s | 成功、`truncated`。ビンの一次線（広いビンと \(2^{-10}\) eV のビン） |
| A2 | K_A2、\(N=4\)、D1000、E_A1 | 成功、`truncated`、`support_order` 30。行内: 返ったパネルのラプラス変換（\(\tau=0.05\)、0.3 eV⁻¹）と 0 次・1 次のモーメントを閉形式と比べる。\(h_{n,j}\) を閉形式と比べる |
| A3a | K_A3、\(N=4\)、D1000、E_A3 | 成功、`truncated`。実際に負の `hole` のビンがある |
| A3b、A3c | A3a で \(Q_h=0\)、2 | 成功。A3c は `negative_total_kernel` の警告（\([1.5,2)\) で `total` が負） |
| A3d | A3a で `multiplier="external"`、\(m=3\) | 成功。\(n\) 次は \(m^n\)。`multiplier_differs_from_weight_length` の警告 |
| A3e | A3a で `multiplier="inverse_window_area"` | 成功。\(w=1\)、\(r=12/19\) |
| A3f | K_A3、\(N=3\)、面積だけ、S2、`inject` で**補間点だけ**に全 11 列の点の誤差 \(10^{-3}\)（値はそのまま、検証点は元のまま）、緩い設定 | 成功。1次: \(A_{\rm total}=0.19\)、誤差は \([0.007,\,0.007(1+10^{-9})]\) Å⁻¹。乗数の区間の両端は、報告した \(A\)、\(e\) から製品の外で作り直した \(1/(A+e)\)（下へ）、\(1/(A-e)\)（上へ）とビット単位で同じ。次数 1 の面積は §4.2 の写像（`total` の返す誤差は旧版の \(|\hat w_1-1|\) で、内訳は `convention` だけ。再帰の種 `fixed_multiplier_area_error` は別に \(\hat me_A\) と照合する）。2・3 次の返す面積の誤差が \(e_f+\rho_n(|w|+e_f)\)。固定の乗数の継承の分は、2次の `total` で \(\hat m^2(0.12\times0.007+0.21\times0.007+0.007^2)\)、下回らず、上は相対 \(10^{-9}\) 以内。3次の `total` の3つの内訳には独立な参照を置かず、作り直しで照合する: 伝播の分は、同じ単位に直した2次の局所の分と \(\|R\|_1\) から。継承の分は、Young の解析の項（\(\|R\|_1e^{(i)}_2+\|\hat C_2\|_1e_R+e^{(i)}_2e_R\)）に、2次の局所の分と \(e_R\) の混合の項を1回だけ足したもの。局所の分は `mesh.local_area_terms` の和。行内の対照: 同じ \(\hat m\) を `external` で与えた実行と、2・3 次の固定の乗数の内訳が等しい（相対 \(4\epsilon\)）。行内の対照（点の経路、X13）: K_A1s、`initial_edges=[2, 5]`、補間点だけに 13a・`no_hole`・`total` の点の誤差 \(e=10^{-3}\)、同じ緩い設定、D1000、E_A1s、\(N=2\)。`_fold_order` に渡った節の内訳を読み、必須の節 \(s=7\) eV での `total` の継承の分（\(m\) を掛ける前）が \(3(0.1e+0.1e+e^2)=0.000603\) と絶対の幅 \(10^{-12}\) で合う |
| A3z | K_A3z、\(N=3\)、D1000、E_A3 | 成功。`hole` の代数の面積は 0（床の内）でビンは非ゼロ。行内: `hole` の列の 2・3 次のモーメント |
| AE1 | K_LIN、\(N=2\)、面積だけ、`initial_edges=[2, 4, 5]`、`inject` で補間点だけに 13a・`no_hole`・`total` の点の誤差 \(10^{-3}\)、緩い設定 | 成功。2次の `hole` の継承の面積の誤差が \(E=m^2\times0.05\times0.003=6\times10^{-4}\) を下回らず、上は相対 \(10^{-9}\) と明示した丸めの内（\(|{\tt hole}|\) を零点 \(T=3\) で分ける）。局所の分は、`mesh.local_area_terms` の和からの作り直し（1 キー）と、局所の予算以下であること（1 キー）を別々に照合する。前提: 1次のパネルが \([2,4]\)、\([4,5]\)（違えば `not_run`） |
| OI1 | K_A1s、\(N=2\)、D1000、E_A1s、`order_inject` で2次の `total`・`no_hole`・13a の補間点に局所の点の誤差（幅 3 eV 以上のパネルで \(10^{-3}\)、幅 1.5 eV 以下で \(10^{-12}\)。パネルの単位） | 成功。初めの網 \([4,7]\)、\([7,10]\) は (iii) を超える（`mesh.history`）。二分が 1 段以上で、最終のパネルはすべて幅 1.5 eV 以下。`bin[1]`（KE 990–994、損失 6–10 eV。\(s=7\) を含む）で、局所の分が正、伝播の分が厳密に 0、継承の分は \(10^{-10}|y|+10^{-13}\) 以下 |
| OI2 | K_A1s、\(N=3\)、D1000、E_A1s、`order_inject` で2次の全パネルの同じ3列に局所の点の誤差 \(10^{-12}\) | 成功、`truncated`（\(N=3<N_D=6\)）。4つの不等式は、2・3 次が届く `bin[1]` の `total` で見る（注入しない `hole` には何も求めない）: 3次の伝播の分 (ii) が正。3次の継承の分 (i) は `A1s.result` と相対 \(10^{-6}\)＋\(10^{-18}\) で合う（注入した局所の分との微小な混合の項が入る）。2次の局所の分が正、2次の伝播の分が厳密に 0。2次の最終の網は `A1s.result` と同じで、`mesh.builds` は 1（2次に入るのは1回。前の網に戻らない） |
| A4a、A4c | 層C。`Panel.from_lobatto(2, 5, …)` の1枚（K_A4_9、K_A4_16 の節の値、誤差 0）、\(N=2\) | 成功。2次のパネルは \(s=7\) に境を持ち、そこの端の節の値（13a、`total`）が \(0.05^2\times3\times18/19\)、\(0.05^2\times3\times34/33\)（相対 \(10^{-13}\)）。9 点則なら前者は \(0.05^2\times3\) |
| A4b | K_A4_9、\(N=3\)、D1000、E_A1s | 成功。2次は 19 次の区分多項式で、補間し直しは厳密でない。2次の二分が 1 段以上 |
| A5 | K_A5、\(N=3\)、D1000、E_A5 | 成功。2・3 次の二分が 1 段以上 |
| I-shift | A1s の一次線と端を +8 eV（D1008、E_A1s_plus8） | 成功。全次数の値と誤差が `A1s.result` とビット単位で同じ |
| I-areas | K_A3、\(N=4\)、面積だけ | 成功。\(w_n\) が `A3a.result` と2つの誤差の和の内で合い、厳密値とも合う |
| INV1、INV2 | K_A3、\(N=2\)、5、D1000、E_A3 | 成功。合計 `spectrum` を、参照の生成器が逆変換の式（§6.1）から作った値と比べる。末項は損失 \(\ge1+2N\) eV（INV1 は 5 eV、INV2 は 11 eV）に届く。境より完全に手前のビンは正常の対照で、一致だけを求める（末項は厳密に 0）。境をまたぐビンは届く側に入れる。届く側では、参照の担当が、末項を外した値が幅を超えるビンを 1 つ以上記録する。両側の外を含む全領域で製品の合計と照合する |
| CT1 | K_A3、\(N=3\)、D1000zero、E_A3 | 成功、`complete_on_display`、`support_order` 0。ビン・外は厳密に 0、面積は計算する |
| CT7 | K_A2 を窓 \([1.3(1+2^{-8}),7.7]\) で、\(N=3\)、面積だけ | 成功。行内: `A2.result` の \(h_2\)、\(h_3\) が \((w_n(a)-w_n(a'))/\ln(a'/a)\) と 2 % 以内 |

**参照の模型と打ち切り**

| ID | 入力 | 期待 |
| --- | --- | --- |
| REF1、REF2、REF3 | 参照の入口。K_A3、\(N=3\)、D1000、E_A3。模型は幾何、階乗、C | 成功、`reference_only`。REF2 は REF1 の \(1/n!\)。REF3 の2次と `A3a.result` の2次の差が \(\tfrac12\bar g^{*2}\)（`total`、`no_hole`） |
| REF4a | 解析の入口（B1）。K_R0、\(N=1\)、D1000、E_A1s | 成功、`truncated`、`support_order` 6。全ビンの上界と `neglected_total_area` が厳密に 0 |
| REF4b、REF4c | 参照の入口。同じ入力で幾何、階乗 | 成功。模型の上界が、\(K_1*K_1\) の届くビン（損失 4–10 eV）で正で、その模型の厳密な落とした量以上 |
| T1a、T1b、T1c | K_T1、\(N=2\)、D1000、E_T1_eq、E_T1_up、E_T1_down | T1a: `truncated`、`support_order` 3（等号の次数を残す）。T1b: `complete_on_display`、`support_order` 2、上界は厳密に 0。T1c: `truncated`、3 |
| T2a | K_T2、\(N=4\)、D1000、E_T2、許容 0.05 | 成功、`bounded_on_display`、`support_order` 50 |
| T2b | 同じ入力で許容 \(2.5\times10^{-3}\) | `truncation_bound_exceeds_tolerance`、`required_order_estimate` \(\ge5\)。許容は `bin[1]`（\(U=30\) eV）と `bin[0]`（\(U=50\) eV）の上界の間にある |
| T2c | 同じ入力で許容なし | 成功、`truncated`。`bin[0..3]` の上界を、参照の「格子の最小」と比べる（下側 \(1-10^{-9}\) と上側 \(1+10^{-6}\) を別々に判定）。上界 \(\ge\) 厳密な累積の落とした量（別の判定）、`bin[0]`・`bin[1]` で上界 \(\le10\times\) 落とした量。`bin[0..3]` は \(\tau>0\)。`bin[4]`・`bin[5]`（\(U_D<5a\)）は `(0.0, "complete", None)` で、\(\tau\) は求めない。`neglected_total_area` の区間が \(Swr^4/(1-r)\) を含む |
| T2e | T2c に `inject` で**補間点だけ**に点の誤差（13a \(1.44\times10^{-4}\)、15-T1 \(10^{-4}\)、`no_hole`・`total` \(2.44\times10^{-4}\)）、緩い設定 | 成功。`bin[0..3]` の上界が、誤差の項 \(e^{-\tau t_{\rm lo}}e_{\rm パネル}\) を入れた参照の値（最適の \(\tau\) を選び直す）と同じ幅で合う。参照の担当は、誤差の項のない値との差が幅を超えるビンを記録する。前提: 1次のパネルが 1、2、4、8、16、32、50（違えば `not_run`） |
| T4u | K_T4_one、\(N=3\)、D1000、E_A1s、`inject`（両方の役割）で 13a・`no_hole`・`total` に点の誤差 \(10^{-9}\) | 成功。\(\hat r\) の区間が 1 をまたぐ。`neglected_total_area` は `None`、理由 `undetermined`。\(U_D=12\)、10、6、5、4、2、0 eV、\((N+1)a=4\) eV なので、`bin[0..4]` は有限の上界と \(\tau>0\)（等号の `bin[4]` は届く側）、`bin[5]`・`bin[6]` は `(0.0, "complete", None)` |
| T4d、T4v | K_T4_above、K_T4_below | 成功。T4d は `None`、`diverges`。T4v は値を返し、区間が \(1024(1-2^{-10})^4\) を含む |
| T5、T5n | K_T5、\(N=2\)、D1000、E_T2、`multiplier="external"`、\(m=10^{25}\)。T5 は許容 \(10^{-3}\)、T5n は許容なし | T5: `truncation_bound_unavailable`。T5n: 成功、`truncated`、`bin[0..3]` は `(None, "no_usable_tau", None)`、`bin[4]`・`bin[5]`（\(U_D<3a\)）は `(0.0, "complete", None)`、`neglected_total_area` は `diverges` |
| DS0 | 公開。`synthetic_wide`、`"display_support"`、D998、E_DS0（一次線が表示より低い KE にある） | 成功、`complete_on_display`、`support_order` 0、`computed` 0、`orders` は空のタプル、M4a を1回実行して配送 1 回（どちらも終了状態の必須の属性）。全領域で `spectrum` は `zero_loss` と同じ（\(J_0=1\) は低い側の外）。表示の上界は 0。`neglected_total_area` は \(N=0\) の式（\(Sw/(1-r)\) と区間。\(w\) の区間は \([m_-,m_+]\) と \([A_{\rm total}-e,A_{\rm total}+e]\) の四隅、\(\hat r\) の区間は §4.2 の規則で作り、\(Sw/(1-r)\) を両端の隅で厳密な有理数で評価して、最後に外向きに丸める） |
| DS1 | 解析。K_T1、`"display_support"`、D1000、E_T1_eq | 成功、`complete_on_display`、`computed` 3、`support_order` 3（等号の次数を計算する）。等号で残した3次は表示と接するだけで、正の質量は低い側の外にある。3次数のビン・両側の外（54 キー）と面積（33 キー）を丸めの水準で比べる。上界は全ビン `(0.0, "complete", None)` |

T2 の許容は、下書きの担当が厳密な \(\Lambda\) で求めた上界（`bin[0]` 約 \(1.02\times10^{-2}\)、`bin[1]` 約 \(5.2\times10^{-4}\)。`work/consult/m4b/prereg/t2_bound_check.py`）を挟むように選んだ。この数は期待値ではなく、参照の担当が R-A で作り直す。

**実材料（公開の入口。`epsrel` \(10^{-7}\)、\(N=3\)、E_B120）**

| ID | 条件（共通の条件から変えた点） | 期待 |
| --- | --- | --- |
| RB1、RB2 | Al、D1000、`k_total`・`k_parallel` | 成功、`truncated`。RB2 は `dispersive_k_parallel_13a` |
| RB3 | Al、B5 | 成功 |
| RB4、RB5、RB6 | Si（窓 \([5,50]\)）、D1000 の `k_total`・`k_parallel`、B5 | 成功。M4a の `resolution_basis_computed_candidates_only` の警告を継承（必須の警告） |
| RB7 | Al、D1000、`multiplier="external"`、\(m=40\) Å | 成功。`neglected_total_area` は `None`、`diverges`。`order_area_ratio_not_below_one` |
| RB8 | Al、D1000、\(Q_h=0\) | 成功 |
| RB1a | Al、面積だけ | 成功。\(w_n\)、\(h_n\)、\(\hat r\) が `RB1.result` と2つの誤差の和の内で合う（網が違いうるのでビット一致は求めない）。1次と `single_loss` は同じプロセスの `loss_window_area` とビット単位で同じ。`mesh.cells` と `formal_folds` は 0 |

全行で、1次と `single_loss` が同じプロセスの `single_loss_spectrum` とビット単位で同じことと、`_call_kernel` の実引数を確かめる。2・3 次のビン、両側の外、面積を R-B と比べる。

**M4A（基準のコードとのビット一致。`synthetic_wide`、窓 \([1,50]\)、D1000、E_M）**

| ID | 規約（M4a-2 設計 §2.4 の行） | 期待 |
| --- | --- | --- |
| M4A-c1、c2 | 旧経路の2行: 省略／`external` で \(m=20\) | \(N=1\)。`single_loss` の全欄（旧来の欄、`diagnostics`、警告）と次数 1 が、基準のコードの `single_loss_spectrum` とビット単位で同じ。`_call_kernel` の呼び出しが同じ。今のコードの `on_final_panels` の配送が 1 回（基準の子には配送の口がないので、今のコードの側で `_run` を包んで見る） |
| M4A-c3–c7 | 追加の5行: `external` 30／S2／窓の 13a の重み（\(c=0.01\) Å⁻¹）と `weight_length`／同と `external` 30／同と S2 | 同じ |
| M4A-o3、M4A-o3s | c1、c4 の規約で \(N=3\) | 同じ（1次の欄は \(N=1\) の行と基準の両方とビット単位で同じ）。高次の成分の恒等式。各高次の正式の畳み込みは 1 回以上。`on_final_panels` が受け取ったパネルのオブジェクトを直接に捕まえ、それが2次の1次の側の入力（`panels_by_order[1]`）であることを見る |
| M4A-Al | Al、c1、\(N=1\)（M4a の RB1 の条件） | 同じ（全件だけ） |
| M4A-Z | AlZ、c1、\(N=2\) | 成功。評価 0 回。高次の値・誤差・面積は厳密に 0、`spectrum` は `zero_loss` と同じ |
| M4A-area-c1、-c4 | 面積だけの入口、\(N=1\) | 旧来の面積の欄と次数 1 の面積が、基準の `loss_window_area` とビット単位で同じ |

- **基準のコード:** 公開したタグ `v0.3.0`（コミット `c40a38940d52f5e637ee987b9c6aad911adca4dc`）。設計 §10 の `5a9f9b3` とは `src/` の中で `__init__.py` の版の文字列の1行だけが違う（`git diff --stat 5a9f9b3 v0.3.0 -- src` で確かめた。Codex も確かめた）。タグは動かず、公開した配布物と同じなので、こちらを採る。凍結の記録に、設計の旧いコミットとの関係を1行書く。
- **方法**（M4a-2 の V1 と同じ）: `git rev-parse v0.3.0^{commit}` を上の値と照合し、コミットを固定した `git archive` で `src/xpsloss`（同梱の比較の表を含む）を取り出す。基準の側は、新しい API を読み込まない小さな入力の変換器で、M4A の全行をまとめた1つの子プロセスで動かす（M4A-o3・o3s の同じ旧版の入力は再利用する）。`sys.executable`、NumPy、SciPy、数値の環境を合わせ、`xpsloss.__file__` と呼び出すモジュールが展開先から来たことを確かめる。ソースの試験でタグや履歴がなければ失敗。展開した sdist の既知の skip は比較の件数に入れない。
- **比べ方:** 構造体の全部の欄の名前と型、辞書のキー、警告の順、配列の dtype・形・バイト列、スカラーの binary64 のビット（符号つきのゼロを含む）、例外の公開の属性を再帰して比べる。旧版の面積（Å⁻¹）と新しい無次元の面積を、同じ裸の数として比べない。`single_loss` は旧版の結果の全欄と照合し、次数 1 の派生の欄には §4.2 の写像を使う。今のコードの各入口は1回だけ実行し、配送を確かめる経路をキャッシュで飛ばさない。
- **版の欄:** `run.json` の `package_version` と `capabilities` の `version` だけは、値で比べず、各プロセスが自分の `__version__` を報告することを照合する（M4a-2 の凍結の後の変更の教訓）。基準の側は、タグのソースの値 `0.3.0` とも独立に照合する。`schema_version`、模型の ID、依存のライブラリの版は、欄ごとの本来の検査をする。基準のコードに新しい CLI はないので、C01 は今のコードの自己報告、C11 は両側の自己報告を別々に確かめる。時刻、所要時間、一時のディレクトリは妥当性だけを別に見る。
- 既存の試験と参照（M4a の 100 行、M4a-2 の行）は変えずに合格すること（行は足さない）。

**拒否・未収束（公開の行は `synthetic_wide`、窓 \([1,50]\)、D1000、E_M、\(N=3\)）**

| ID | 入口と入力 | 期待 |
| --- | --- | --- |
| N01、N02 | `model` が `"full_kernel_geometric_ref"`／`"b1"` | `invalid_model`、評価 0 回 |
| N03–N06 | `max_order` が 0／2.5／`True`／`"auto"` | `invalid_max_order`、評価 0 回 |
| N07 | 面積だけの入口で `"display_support"` | `invalid_max_order` |
| N08–N10 | 許容が −1／`nan`／0 | `invalid_neglect_tolerance` |
| N11 | `max_order` 65 | `order_limit_exceeded`（`requested_order` 65、`order_limit` 64） |
| N12 | `"display_support"`、E_D100（\(N_D=100\)） | `order_limit_exceeded`（`requested_order` 100） |
| N13 | \(N=4\)、`max_convolution_points` 4500 | `work_lower_bound_exceeded`（`resource` は `max_convolution_points`、`lower_bound` 4743、`limit` 4500） |
| N14 | \(N=40\)、`max_cells` 300 | `work_lower_bound_exceeded`（`max_cells`、`order` 34、`lower_bound` 306）。M4a のセルの検査より前 |
| N15 | \(N=6\)、`max_convolution_points` 50,000 | `convolution_work_limit_exceeded`、`stage` は `after_first_order`、評価 1 回以上。入力の下限 10,540 は通り、見積りは \(P_1\ge9\) で 70,125 以上。行内: 報告した `estimate` が式と報告した \(P_1\) に合う |
| N16 | 解析。K_STEP、\(N=2\)、D1000、E_A1s、`max_convolution_points` 4500 | `convolution_work_limit_reached`、`stage` は `order_refinement`、`order` 2。入力の下限 1054 と、\(P_1=6\) の見積り 3927 は通る。2次の必須の 12 枚の内部の点だけで \(12\times31\times17=6324\) 点が要る。止まった時点の実測（`measured`）を返す。前提: 1次が 6 枚（違えば `not_run`） |
| N18 | 解析。K_A5、\(N=2\)、D1000、E_A5、`order_panel_limit` 4 | `order_budget_not_met`、`order` 2（枚数の上限で止まる） |
| N19 | 解析。K_SGN、\(N=2\)、D1000、E_SGN_rel、`inject`（両方の役割）で 13a・`no_hole`・`total` に点の誤差 \(1.5\times10^{-8}\)、`spectrum_epsabs_fraction` \(10^{-300}\)、`area_epsabs` \(10^{-6}\) | `inherited_error_exceeds_order_budget`、`order` 2。`exceeded` に `bin[3].hole` があり、その `side` は `relative`、`suggestion` は `epsrel`。面積のキーは `exceeded` にない。`bin[3]`（損失 6.875–7 eV）で、継承の包絡の積分が \(|y|\) の \(7.05\times10^{-6}\) 倍（予算は \(2\times10^{-6}\)）。その要素の `inherited`、2つの予算、採った側を照合する |
| N20 | 解析。同じ核、E_SGN_abs、点の誤差 \(10^{-8}\)、`area_epsabs` \(10^{-6}\)、ほかは既定 | 同じ理由。`exceeded` に `bin[2].hole` があり、その `side` は `absolute`、`suggestion` は `spectrum_epsabs_fraction`。面積のキーは `exceeded` にない。`bin[2]`（損失 6.5–7.5 eV）は反対称で \(y=0\)、継承の包絡の積分 \(4.4\times10^{-9}\) が絶対の予算 \(\phi_2\eta S_{\rm abs}=8.4\times10^{-10}\) を超える |
| N21 | 解析。K_A5、\(N=2\)、`order_min_panel_width_eV` 4.0 | `order_mesh_unresolved`、`order` 2、`panel_eV` つき |
| N26 | 解析。K_A1s、\(N=2\)、`external`、\(m=10^{200}\) | `nonfinite_output`、`order` 2（1次は有限） |
| N27 | 窓 \([0,50]\) | M4a の `SingleLossInputError`、`invalid_window` がそのまま通る |
| N28、N29 | AlZ と S2／AlZ、窓の 13a の重み（\(c=0\)）と S2 | `multiplier_undefined`／`weight_length_undefined`（先に出る） |
| N30 | `panel_limit` 9 | M4a の `SingleLossConvergenceError`、`loss_budget_not_met`（`dominant_part` は `interpolation`）。高次は始まらない |
| N31、N32、N33 | `order_epsrel` が `nan`／`max_convolution_points` が `True`／`order_limit` が 0 | `invalid_settings`（`argument` つき）、評価 0 回。表の `nan` の印は実際の浮動小数に戻して渡す。真偽の型は整数への変換の前に検査する |

N19 と N20 の1次は通る（Codex の実測で、1次の予算に対する比の最大は 0.50 と 0.33）。`area_epsabs` を \(10^{-6}\) にしたのは、既定のままだと2次の `hole` の面積の絶対の予算（\(1.68\times10^{-10}\)）を Young の項（\(4.3\times10^{-8}\)、\(2.9\times10^{-8}\)）が先に超え、意図したビンの側が選ばれないためである。\(10^{-6}\) なら面積の予算は \(1.68\times10^{-6}\) で、面積は通る。継承の見積りは、包絡の積分（上の数）以上になる（セルの重みは補間の基底の絶対値なので）。試験は代表のキーの順に依らず、`exceeded` の中の対象のビンを見る。

**CLI（`synthetic_wide`、\(E_0=1000\) eV、窓 \([1,50]\)、\(L_w=20\) Å、`k_total`、`--model first_loss_13a_chain_v0 --max-order 3`、D1000、E_M）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| C01 | 上のまま | 完了、終了 0、`truncated`。行内: `effective_settings`（13 個）からの API の再生と、全次数の値・誤差・2つの CSV の全数値の列がビット単位で一致。`run.json` の記録の欄。同じ入力の再生で入力のハッシュ、実効の設定、状態、警告を照合 |
| C02 | `--max-order display_support`、E_DS2 | 完了、`complete_on_display`、`computed` 2、`support_order` 2 |
| C03、C04 | `--model`／`--max-order` を与えない | `error`、`ValueError`、`reason` なし、終了 2、出力先なし |
| C05 | `--max-order 1 --neglect-tolerance-fraction 1e-30` | `not_converged`、`truncation_bound_exceeds_tolerance`、終了 2、出力先なし。stderr に `order`、`stage`、`inner`、`required_order_estimate` |
| C06 | `--max-order 1 --multiplier-A 1e27 --neglect-tolerance-fraction 1e-3` | `not_converged`、`truncation_bound_unavailable`（C05 と別に照合する） |
| C07 | 既存の出力先 | `error`、`FileExistsError`、終了 2。既存のファイルのバイトは変わらない |
| C08 | `--max-order 65` | `error`、`order_limit_exceeded` |
| C09 | `spectrum.json` の `axis` が `"binding_energy_eV"` | `error`、`MultipleLossInputError`、`invalid_schema` |
| C10 | `--max-order 6 --max-convolution-points 50000` | `not_converged`、`convolution_work_limit_exceeded`、stderr の `stage` は `after_first_order`、`estimate` と `limit` つき。生成物なし |
| C11 | `xpsloss capabilities` | 終了 0。§4.3 の `capabilities` の項目 |

**対照（公開の入口）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| CT3 | `synthetic_wide`（無分散）、`k_parallel`、\(N=3\) | 成功。`dispersive_k_parallel_13a` は出ない。全次数の値が `M4A-o3.result`（`k_total`）と2つの誤差の和の内で合う |
| CT4、CT5、CT6 | Al、E_M、\(N=1\)、窓の 13a の重み、`weight_length`、\(c=0\)、\(10^{-2}\)、\(10^{-10}\) Å⁻¹（全件だけ） | 成功。`convention_ratio` と `neglected_total_area` が `exactly_one`＋`undetermined`／`below_one`＋値／`below_one`＋`undetermined` |

設計 §10 の対照のうち、\(K_1=R\) は A1、\(R=0\) は REF4a、\(F=0\) は CT1、構造的なゼロは M4A-Z・N28・N29、\(m\ne L_w\) は A3d・RB7、\(h_n\) と窓の端の差分は CT7 が持つ。

### 5.3 CI と全件

- **CI（90 行）**と**全件だけ（19 行）:** RB2–RB8、RB1a、A5、CT3、M4A-c2、-c3、-c6、-c7、-Al、-area-c4、CT4–CT6。A5 と CT3 は初めから全件だけにした（Codex の提案）。
- CLI の行は同じプロセスの中で `main()` を呼ぶ。基準のコードの側は1つの子プロセスで行う。
- **CI の増分の測り方:** 製品の実行、同じプロセスの中の M4a の対照、基準の子プロセスの起動をすべて含めた経過時間で測る。試走の CPU 秒と、並行して動かしたときの経過の秒を足し合わせない。
- **60 秒を超えたときの順（表の `ci_reduction_order`）:** M4A-o3s、C02、RB1 の順に1行ずつ全件だけに回し、回すたびに測り直す。C01（独立な API の再生）、N13・N15（段の区別）、M4A-c1・c4・c5・o3 は残す。C02 を回しても、`"display_support"` の成功の軽い行（DS0、DS1）が CI に残る。回した行は全件と変異の実行の対象に残し、各段階の件数と再利用の向きを生成で検査し直す。それでも超えれば、旧版を残して再登録する。
- 全件の記録は、統括担当が全 109 行を `execution_order` で1回実行し、終了状態、数値の精度、推定の被覆、参照の経路の失敗、キーの欠落・重複、未実行、変異を別々の欄で検証記録に書く。

## 6. 参照値（実装担当とは別の担当）

参照の担当は、実装のコードと結果を見ずに作る。実装担当は参照値を見ない。参照値は凍結して SHA-256 を付録に記録してから、実装と比べる。値は十進 30 桁以上の文字列で保存し、\(u_r\) に十進への丸めの幅を含める。**恒等式の両辺を製品の畳み込みで作らない。**

### 6.1 R-A（解析の参照）

- **パッケージ `xpsloss` を読み込まない。**二進の入力を厳密な有理数として読み、有理数か mpmath（40 桁以上）で評価する。
- **箱形の核（K_A1、K_A1s、K_A3、K_A3z、K_T1、K_T2、K_T4、K_R0、K_SGN）:** 幅の違う一様変数の和の累積分布（分母 \(n!\prod h_i\)、端のずれの部分集合の和、符号は要素数の偶奇）を有理数で厳密に。δ の一次線のビンは累積分布の差、ビンの一次線はその積分（有理数）。両側の外は、その領域に入る部分を直接に求める。
- **K_A2:** \(C_n(s)=\kappa^ne^{-\mu s}\times\)（ずらした Irwin–Hall の密度）。ビンは多項式×指数の積分の閉形式。ラプラス変換とモーメントも閉形式。\(h_{n,j}\) は \(h_1r^{n-1}+(n-1)w\,h_rr^{n-2}\)。
- **K_A4:** 多項式の畳み込みを有理数で厳密に。\(s=a+b\) の点の値は \(\kappa^2\ell\cdot18/19\)、\(\kappa^2\ell\cdot34/33\)。
- **K_A5:** 多倍長の入れ子の求積（\(n\le3\)）。作業精度と次数を変えた差を \(u_r\) に入れる。
- **INV（R-INV。Codex の案で確定）:** 損失の座標で累積 \(H(s)=\int_{-\infty}^{s}J^{(N)}\) を解く。\(H(s<0)=0\)、原点に δ の質量 \(S\)、\((I-\mathcal A_R)J^{(N)}=(I+\mathcal A_{\bar g})F-\mathcal A_{K_1*R^{*N}}F\) の \(\mathcal A_R\) は箱形の \(R\) による遅延の積分になる。遅延が 2 eV 以上なので、既知の側から幅 2 eV ずつ区分的に閉じる。右辺の打ち切りの項 \(K_1*R^{*N}\) は箱の独立な式で作る。製品の `spectrum` のビンと両側の外を、この独立な解の累積の差と台の全体の積分と比べる。次数の和による自己検査は、参照の値を作る関数と別の関数にする。末項の台は損失 \(1+2N\) eV から始まり、遅延の積分も差を損失の小さい側に戻さないので、境より完全に手前のビンでは末項を外しても値は厳密に同じである。そこは正常の対照（一致だけ）とし、届く側（境をまたぐビンを含む）で、末項を外した誤った値が幅の外にあるビンを 1 つ以上記録する（INV1 は 5 eV、INV2 は 11 eV が境）。幅のある核は製品のビンの移動だけでは表せないので、設計 §10 の「製品の出力に \(\mathcal A_R\) を作用させる」形は取らない。
- **結合則:** \(C_3\) を \((K_1*R)*R\) と \(K_1*(R*R)\) の両方で作り、有理数で等しいことを記録する。
- **上界（T2c、T2e、REF4b、REF4c）:** §3.3 の式を登録した格子で評価した最小（T2e は1次のパネルの端 1、2、4、8、16、32、50 と注入した誤差から）。厳密な落とした量は Irwin–Hall の有限和。
- **誤差の内訳（A3f、AE1）:** Young の右辺を厳密な \(L^1\) のノルムで作る（3次の \(\|\hat C_2\|_1\) は区分1次式の絶対値の積分）。
- **既知の誤答:** §9 の変異に当たる誤った値を参照の側で作り、§7.2 の幅の外にあることを記録する（X01、X02、X04–X10、X12、X14–X16、X18–X20、X23、X33）。

### 6.2 R-B（実材料）

- **方式（設計 §10、§14 の T11）:** 公開の `effective_cross_section_depth_averaged`（精度は API の既定 \(10^{-12}\)、\(10^{-9}\)、200）の点の値を、窓の端を格子点に置いた一様な格子（\(h=2^{-k}\) eV。\(a\)、\(b\)、ビンの端が格子点に乗る）で取り、台形則で畳み込み、\(h\) を半分ずつにして Richardson 外挿する。端の和 \(2a\)、\(a+b\)、\(2b\)（3次は 4 点）は別に扱う。ビンは外挿した累積分布の差（ビンの一次線はその積分）。面積は窓の全体で別に積分する。
- **独立性:** 公開の点の値だけを共有する（共有する部分は照合の範囲の外と明記する）。実装のパネル、畳み込み、セルを使わない。`total` は M3a の `total` を直接に使う。
- **\(u_r\) の内訳**（キーごとに記録して足す）: (1) 最も細かい2つの \(h\) の外挿の差、(2) 外挿の回数を変えた差、(3) 上流の点の誤差を畳み込みの重みの絶対値で伝えたもの、(4) ビンの一次線の積分の差、(5) 丸め。\(r\ne0\) のとき \(f_r=u_r/|r|\)。
- **比較の前に記録すること:** 2つの運動量の指定の参照の差が、2次の主なビンで比較幅の和を超えること（X01 の誤った値が RB1 の幅を超えることも）。

### 6.3 「値なし」の基準（実装の値を見る前に凍結する）

1. **\(r=0\) のキー:** 床 \(d_{\rm abs}\) が正なら \(u_r\le d_{\rm abs}/4\) を確かめて値を持つ。厳密な 0（`exact_zero_or_null`）は \(u_r=0\) と \(y=0\) を求める。
2. **\(r\ne0\) の R-A:** \(u_r\le\)（§7.2 の幅）\(/4\) を満たさないキーは `no_reference`。
3. **\(r\ne0\) の R-B:** \(f_r\le1.25\times10^{-5}\) を超えるキーは `no_reference`。
4. 参照の経路が失敗したキーは「参照の経路の失敗」として別に数える。
5. **必須のキーが得られない行を成功にしない:** 行の外部の参照のキーの 1 割を超えて `no_reference` か経路の失敗なら、その行は「参照不足」と記録し、合格に数えない。`direct_reference` が 0 件の行（AE1、T2e、REF4b、REF4c など）は、唯一の必須の参照（`error_part` や `bound_reference`）が1つでも欠ければ参照不足。RM-B の**スペクトルの行**は、識別のキー（2次の `total` の最大のビンと、その両隣）が1つでも欠ければ参照不足。面積だけの RB1a は、2・3 次の `total` の面積と `total` の感度が必須。

`no_reference` のキーは合格にも不合格にも数えず、台帳に残す。「照合未確定」は M4a §6.3 の手順（キーと必要な精度だけを参照の担当に渡す）で扱う。**別々に数えるもの:** 合格、不合格、`no_reference`、`unresolved`、`not_run`、実装の予期しない未収束、上流の失敗、参照の経路の失敗、参照不足の行。

## 7. 比較のキーと許容幅

### 7.1 キー

- **スペクトル:** `order[n].<領域>.<total|no_hole|hole>`（領域は `outside.low`、`bin[i]`、`outside.high`）、`spectrum.<領域>`。
- **面積:** `order[n].area.<列>`、`order[n].cutoff_sensitivity.<列>`、`order_area_ratio`。
- **誤差の内訳:** `order[n].<成分>_error_parts[part][i]`、`order[n].outside[side]["<成分>_error_parts"][part]`、`order[n].area_error_parts[part][column]`（`part` は `inherited`、`propagated`、`local`、`multiplier_interval`。面積の内訳だけ `convention` を加えた5つ）。
- **打ち切り:** `truncation.neglected_at_or_above_region_bound[bin[i]]`、`truncation.neglected_total_area`、`truncation.state`、`support_order`。
- 合否は「行、キー」の単位で数える。`total`、`no_hole`、`hole` は別に判定する。部分的に重なるビン（参照と端が違う）は置かない（参照は行の端そのもので作る）。登録した参照のキーの未使用と重複の使用は失敗にする。

### 7.2 精度の比較（M4a §7.2 の形を次数に広げる）

\[
|y-r|+u_r\le d_{{\rm abs},n}+d_{{\rm rel},n}|r|
\]

**この幅は、登録した入力での、独立な参照に対する精度の目標である。**実装の返す誤差と独立に決める。伝播の分 (ii) も乗数の区間も、この式で抑えられるとは言えない。だから、API の成功、§7.2 の外れ、§7.3 の過小推定・未確定を別々に数え、返った (ii) を幅に足して通すことはしない。\(r=0\) のキーでは幅を \(d_{{\rm abs},n}\) とし、\(u_r\le d_{{\rm abs},n}/4\) のとき \(|y|+u_r\le d_{{\rm abs},n}\) を検査する。

- **参照の側のノルム:** \(w_{{\rm abs},j}=\|K_{1,j}\|_1\)、\(r_{\rm abs}=\|R\|_1\)（\(m\) を含む無次元。参照の担当が厳密に、または R-B の格子で作る）、\(\phi_{n,j}=r_{\rm abs}^{\,n-1}+(n-1)w_{{\rm abs},j}r_{\rm abs}^{\,n-2}\)、\(\psi_n=\sum_{k=2}^{n}r_{\rm abs}^{\,n-k}\)。\(n=1\) は \(\phi_1=1\)、\(\psi_1=0\) と置き、\(R=0\) で負のべきを評価しない（\(\phi_2=w_{{\rm abs},j}\)、\(n\ge3\) で 0。\(\psi_n=1\)）。
- **標準の幅:**

| 量 | \(d_{{\rm rel},n}\) | \(d_{{\rm abs},n}\) |
| --- | --- | --- |
| ビンと外（\(I\)） | \(n\varepsilon+(n-1)\varepsilon_o\) | \((\phi_{n,j}\eta+\psi_n\eta_o)S_{\rm abs}\) |
| 面積 \(w_{n,j}\)（無次元） | 同じ | \(\phi_{n,j}m\eta_A+\psi_n\eta_o\) |
| 感度 \(h_{n,j}\)（無次元） | 同じ | \(n\,m\,a\times10^{-12}\,\max(1,r_{\rm abs})^{n-1}\)（\(m\) は Å、\(a\) は eV、\(10^{-12}\) は Å⁻¹ eV⁻¹） |
| 合計 `spectrum`（INV。\(I\)） | \(N\varepsilon+(N-1)\varepsilon_o\) | \(\eta S_{\rm abs}+\sum_{n\le N}d_{{\rm abs},n}\)（`total` の列） |

  \(n=1\) では M4a §7.2 の幅になる。行の設定（RB の `epsrel` \(10^{-7}\) など）をそのまま入れる。R-B では、参照が適格なキー（§6.3）だけに、\(d_{{\rm rel},n}\) を \(\max(8f_r,\ n\varepsilon+(n-1)\varepsilon_o)\)、最大 \(10^{-4}\) として当てる。
- **丸めの水準（`exact_polynomial`。A1、A1s、A1p、T2c、DS1 のビンと面積）:** \(d_{\rm rel}=10^{-12}\)、\(d_{\rm abs}=10^{-13}S_{\rm abs}\)（面積は \(10^{-15}\)）。この幅は、ここに挙げた低い次数の箱形の核と、その最終のセルでの回帰の目標に限る。一般の入力の精度の約束ではない。
- **点の値（A4a、A4c。パネルの単位 Å⁻² eV⁻¹）:** 相対 \(10^{-13}\)、床 \(10^{-17}\)。\(s=7\) の両側のパネルの端の節の値を読む。**ラプラス変換と 0 次のモーメント（Å⁻ⁿ）、1 次のモーメント（Å⁻ⁿ eV）:** 相対 \(n\varepsilon+(n-1)\varepsilon_o\)、床は面積の床を \(m^n\) で割った値（1 次のモーメントはその \(nb\) 倍）。これらの派生の量には返す誤差がないので、精度だけを登録し、§7.3 を課さない。
- **上界（T2c、T2e）:** 参照の格子の最小に対する下側 \(1-10^{-9}\) と上側 \(1+10^{-6}\) を別々に判定する。上側は、定数のパネルの指数の積分を安定に求め、左端の重みを誤差の項にだけ掛ける作り方（§3.3）に限って成り立つ目標である。「真の累積の落とした量を覆う」（上界 \(\ge\) 落とした量）は別の判定。`bin[0]`・`bin[1]` の比 \(\le10\) も別の判定。REF4b・REF4c は「正で、落とした量以上」だけ。
- **誤差の内訳（`error_part`。外部の参照）:** A3f の2次の継承の面積の誤差（4 列）、AE1 の2次の `hole` の継承の面積の誤差は、参照の値を下回らず、上は相対 \(10^{-9}\) と、入力の誤差・丸めを明示した分の内。A3f の点の経路の対照は絶対の幅 \(10^{-12}\)。製品が計上する補間の残差と求積の丸め（局所の分）には独立な参照を置かず、作り直し（`recomputed_from_report`）と予算の不等式で見る。
- **区間（`interval_contains`）:** 実装の区間が参照の値（と \(u_r\)）を含めば合格、交わらなければ不合格、ほかは照合未確定。S2 の乗数の区間、\(\hat r\) の区間、DS0 の全軸の残りの区間は、報告した値から §4.2 の丸めの順（四隅を厳密な有理数で求め、最後に外向きに1回丸める）で作り直した両端ともビット単位で照合する（「含む」だけで過大な区間を通さない）。

### 7.3 推定の検査（三値）

返す誤差を持つ `direct_reference` のキーで、\(d=|y-r|\)、\(e=e_{\rm ret}+f_{\rm fp}\)（\(e_{\rm ret}\) は実装が返した誤差、\(f_{\rm fp}=4\epsilon\max(|y|,|r|)\)）として、\(d+u_r\le e\) なら**推定が覆う**、\(d-u_r>e\) なら**過小推定**（不合格）、その間は**照合未確定**。精度の比較（§7.2）とは別の判定で、同じ参照の記録から両方を作って別々に残す。成功した行では、(i) と (iii) の判定が通っていること（`error_parts` の `inherited` と `local` が、返った値から作り直したそれぞれの予算以下）も確かめる。(ii) には予算を課さない。層C、参照の模型、派生の量（ラプラス変換、モーメント、点の値）には、この節を課さない。未使用・重複・空の照合、未確定、未実行は合格に含めない。

### 7.4 行内の恒等式の幅

| 検査 | 幅 |
| --- | --- |
| 成分の関係 | 各次数・各領域・面積で `total` と `no_hole`\(+Q_h\cdot\)`hole` の差が、3つの返す誤差の和と \(8\epsilon\) の丸め以下 |
| 全領域の和 | \(|\sum_Dy_{n,D}-Sw_n|\le\sum_De_{n,D}+e_{w_n}S+8\epsilon(\sum_D|y_{n,D}|+|Sw_n|)\)。合格の根拠にはせず、落ちたら不合格にする |
| 合計 | `spectrum` が \(J_0+\sum_nJ_{n,\rm total}\) と \(8\epsilon\sum|\cdot|\) 以内。`spectrum_error` が返す誤差の和以上で、和の \((1+8\epsilon)\) 倍と \(8\epsilon\sum|\cdot|\) の和以下 |
| RB1a、I-areas、CT3 | 2つの実行の差が、2つの返す誤差の和以下 |
| REF2 | \(|y_{{\rm REF2},n}-y_{{\rm REF1},n}/n!|\le e_{{\rm REF2},n}+e_{{\rm REF1},n}/n!+8\epsilon|y|\) |
| CT7 | \(|h_n-(w_n(a)-w_n(a'))/\ln(a'/a)|\le0.02|h_n|\)（\(n=2,3\)。0 でない列） |
| I-shift、M4A、RM-B の1次、C01 の再生 | ビット単位（配列のバイト列。符号つきのゼロを含む） |

## 8. 時間の上限

- **CI の増分:** 各代表の環境で 60 秒以内（測り方は §5.3）。超えたら §5.3 の順で全件だけに回す。50–60 秒という第1稿の数は見積りだった。Codex の実測では、RB1 の1次だけで、2つのプロセスを同時に動かしたときの経過時間が 16.2 秒である。A5 と CT3 を初めから全件だけにしたので、CI の見積りはその分（約 8 秒）小さい。有理数の端の実装の時間は実装の後に測る。
- **全件:** 30 分以内（試走の見積りは 10–15 分。RM-B の8行で約 130–150 秒）。
- **1行:** 120 秒を超える行は登録しない。一次線は δ か 5 ビン以下、出力は 160 ビン以下に限った。実測の大きさの一次線（試走 T5 の形）は試験の行にしない。
- **参照の生成:** R-A は 5 分以内、R-B は1コアの経過時間で 30 分以内。超えたら比較を開く前に打ち切り、足りないキーを §6.3 で数える。
- **全変異:** 別の枠で、全 29 個の正常版と変異版の対を 60 分以内。超えた変異は「未実行」と記録し、検出に数えない。CI の時間で行を全件に回しても、全件と変異の分母は変えない。

## 9. 変異の表（必須 29）

レビュー担当が、各行の変異を本番のコードに1つずつ入れる。設計 §10 の X01–X24、X28–X31、X33 をすべて必須にした。X25–X27 と X32 は BE 軸の変異なので M4b-2 に回す。対象のキーと枝は表の `mutations` にある。下は実装の前に成り立たせる検査で、検出の実績ではない。

| ID | 変異 | 登録した行と、正常版から区別する観測 |
| --- | --- | --- |
| X01 | 2回目以降の核を `total` にする | A3a の2次の `total`（補助に RB1）。参照の担当が、幅を超える実際のビンの添字を記録する |
| X02 | 2回目以降の核を `no_hole` にする | A3b の2次の `total`（\(Q_h=0\) でも 15-T1 が残るので、13a だけを繰り返す値と分かれる） |
| X03 | 1回目の核を 13a にする | M4A-c1 の旧版との比較、A3a の1・2次のビン |
| X04 | 超過に \(n\) を掛ける | A3a の2次の面積と、非ゼロの `total` のビン（超過が1回分余計に入る） |
| X05 | 各次数に \(1/n!\) を掛ける | A1 の2次のビンと面積（正の2次が半分になる） |
| X06 | 繰り返しに \(m\) でなく \(L_w\) を使う | A3d の2次のビン（繰り返しの因子が 3 から 2 に変わる）、補助に RB7 |
| X07 | 各次数に \(m\) を1回しか掛けない | A3d の2次のビンと面積（\(m^2\) が \(m\) になる）。X06 と別のキーとして登録 |
| X08 | S2 を成分ごとに別々に適用する | A3e の1次の `no_hole`・`hole` の面積と \(\hat r=12/19\)。成分ごとの逆数が例外になるだけの改変は、これらのキーを破らない限り検出に数えない |
| X09 | \(s>b\) の2次をゼロにする | A1s の `order[2].bin[2].total` |
| X10 | \(C_{n-1}(s-t)\) を \(C_{n-1}(s+t)\) にする | A2 の2次のビン、ラプラス変換、1 次のモーメント |
| X11 | 必須の内部の節を1つ落とす（A1 の3次の \({\rm RN}(2a+b)\)） | A1 の3次の、初めの網と最終の網の**実際のパネルの端**（ビット列）。丸めの水準の比較は補助。自己申告の `knots_eV` や数値の精度だけに頼らない |
| X12 | ガウス点を 17 から 9 にする | A4a、A4c の \(s=7\) の両側の節の値（9 点則の検出は A4a が決め手） |
| X13 | 誤差の伝播の項を1つ落とす（面積の経路の3つの枝と、点の経路の3つの枝を別々に実行。ID は1つ） | 面積の経路: A3f の2次の `total` の継承の面積の誤差（3項は \(m^2\) を除いて 0.00084、0.00147、0.000049）。点の経路: A3f の行内の対照で、`_fold_order` に渡った \(s=7\) eV の節の継承の分（0.000603。3つの枝の欠ける量は 0.0003、0.0003、0.000003 で、幅 \(10^{-12}\) より大きい）。点の再帰だけが項を落としても見つかる |
| X14 | A2 の1次を、E_A1 の端を損失の軸に写したセル \([1000-d_{i+1},\,1000-d_i]\) で平均した区分定数（平均の分母はセルの全幅。窓の外は 0）にしてから2次へ渡す | A2 の2次のビン（登録した候補は `bin[5]`）。A2 の両側の外の1次の衛星の質量は 0。参照の担当が、同じ粗視化の既知の誤答で幅の外に出る2次のビンを 1 個以上記録する |
| X15 | 低い側の外を収支の残差で作り、高い側を落とす | A1 の `order[1].outside.high.total`（1012 eV の線から届く）と `order[6].outside.low.total`。両側を直接に照合し、収支だけで判定しない |
| X16 | 全軸の残りを、計算した `outside.low` に足す | A1 の `order[6].outside.low.total` と `spectrum.outside.low`（\(J_0+\sum_nJ_{n,\rm total}\) との恒等式）、T2c の `order[4].outside.low.total`。次数ごとの配列を正しく返して、最後の合計にだけ残りを足す改変は、合計のキーで見つかる |
| X17 | `complete_on_display` の境を狭義の比較にする | T1a（\(N_D=3\)、`truncated`）と T1b（\(N_D=2\)、`complete_on_display`）を一組にする。T1b は正常の対照 |
| X18 | 上界の \(e^{\tau U}\) を落とす／\(\Lambda_R\) の指数を1つずらす（上と下。3つの枝） | T2c の有限で非ゼロの上界と、格子の最小との差 |
| X19 | 和に次数 \(N\) を入れない | A1s の最終の次数が入る `spectrum.bin[0]`・`bin[1]`、INV1 の `spectrum` |
| X20 | 負の値を 0 に切る | A3a の負の `hole` のビン、A3c の負の `total` のビン（符号を含む個別のキー） |
| X21 | 1次を新しい経路で計算し直す／同じ順に2回計算する（2つの枝） | M4A-c1（\(N=1\)）で、旧版との値、`_call_kernel` の呼び出しの列、配送の回数。M4A-o3 で、受け取ったパネルのオブジェクトが2次の入力になったこと（直接に捕まえる）。新しい経路が算術だけを写しても見つかる |
| X22 | \(\rho_n\) の指数を 1 にする | A3f の2・3 次の返す面積の誤差（報告した区間から有理数で作り直す） |
| X23 | \(h_n\) の \((n-1)\) を落とす | A2 の3次の感度（\((n-1)=2\)）。CT7 の差分は別の対照 |
| X24 | 区間を入れた誤差を再帰に戻す（1→2、2→3 の2つの枝） | A3f の2・3 次の固定の乗数の内訳。同じ \(\hat m\) を `external` で与えた行内の対照と等しいことで、広い幅に隠れない |
| X28 | 高次の正式の `total` に \(J_0\) を足す（`_fold` の `spectrum` を読むのと等価） | A1s の `bin[6]`（\(J_0=1\)）での2次の `total` と誤差の厳密な 0、`spectrum.bin[6]=1`。無いキーを読んで例外になるだけの改変は検出に数えない |
| X29 | \(\Lambda_f\) を点の値だけで作る | T2e の `bin[0..3]` の上界。パネルの前提が崩れたら試験は不成立 |
| X30 | 局所の誤差を継承として返し、二分しない | OI1。変異の後は、2次の `inherited_error_exceeds_order_budget`（`bin[1]`–`bin[3]` の `total`）か、二分 0 段の成功のどちらかになる。どちらも、状態・二分の段数・`bin[1]` の内訳を破る |
| X31 | `bounded_on_display` の判定に外の領域を入れる／最も低い KE のビンを落とす（2つの枝） | T2a（外を入れると失敗する。低い側の外に \(\tau\to0\) の極限を使う解釈でも \(1.1956\times0.7056^4/(1-0.7056)>0.05\)）、T2b（`bin[0]` を落とすと成功する） |
| X33 | 絶対値の包絡を、分けないパネルの上でガウス則でそのまま積分する | AE1 の2次の `hole` の継承の面積の誤差（約 \(0.9978E\) になり、\(E\) を下回る）。局所の分を二重に入れた場合は、局所の項の別の照合で見る |

**検出の数え方**（M4a §9、M4a-2 と同じ）:

- 登録した入力で正常版が通り、変異が実際に走り、登録したキー、または予告した型・理由・段の状態の差が出て、初めて1件。**正常版の合格と変異版の不合格をそろって必須にする。**
- 1つの ID に複数の枝があるもの（X13、X18、X21、X24、X31）は、枝ごとにも結果を残し、必須の枝が全部落ちたときに検出済みとする。
- skip、比較 0 件、期待と違う例外、試験器の例外、無関係なクラッシュ、変異が入っていない場合、参照なし、正常版の失敗、登録の外の検査だけの失敗は、検出に数えない。
- 変異の対象の行だけを実行するときも、その行に登録した参照のキーの消費（未使用・重複）を監査する。
- 報告は必須 \(k/29\) とし、生存、未実行、参照不足、正常の対照の失敗、判定不能を別々に示す。BE の4件は次の工程。

## 10. 凍結の前の確認

**下書きの担当が確かめたこと（2026-11-27）:**

- 材料のファイルの SHA-256（Al、Si。記録のファイルと実行用の設定）は `benchmarks/inputs/` の現在のバイト列と一致し、2つのファイルの材料の数値は等しい。
- 基準のタグ `v0.3.0` は `c40a38940d52f5e637ee987b9c6aad911adca4dc` で、`5a9f9b3` との `src/` の差は `__init__.py` の1行だけ。
- 表の生成は、109 行・CI 90・全件だけ 19（登録した定数との照合）、入口・入力・設定の一意性、理由と状態の語彙、公開の理由ごとに1行以上あること（M4b-2 に回した3つを除く。`invalid_settings` を含む）、再利用の向き（先の行だけ、CI から全件だけへの再利用なし）、変異の ID の集合を例外で検査し、§5.1 のキーの群の展開の検査（空・重複、型つきの署名、成果物の供給、変異と台帳の結びつき、CI を回した各段階）を行う。`--check` はバイト単位で一致する。
- 厳密な小さい計算で確かめた数: \(N_D\)（A1 39、A2 30、A1s 6、T2 50、N12 100）、N13・N14・N15・N16 の下限と見積り、N19・N20 の包絡の積分と予算、A3f・AE1 の Young の項、T2 の許容が上界を挟むこと、A1s の台の外の領域の数。
- 現行のコードで確かめたこと: `Panel.check` と `_areas` は、検証点の誤差と補間点の誤差を別の項に足す（Codex の指摘 1 のとおり）。`_fold.propagate_multiplier` は S2 の全軸の重みの誤差を \(|\hat w_1-1|\) にする（指摘 2 のとおり）。解析の核の1次のパネルの数（K_A1 3 枚、K_A4_9 6 枚、K_A5 13 枚）。

**ここで決めた事項（設計が値や名前を書いていない所）:** 既定の高次の予算と資源の上限を API の既定にしたこと。`panel_limit` を高次の網にも当てること。和に予算を課さないこと。求積点の課金の単位と下限の式（§3.4）。\(R\equiv0\) と \(N_D=0\) の状態（§3.1）。ビンごとの上界の3つの場合と、繰り返しの核が負のときの全軸の残り（§3.3）。`neglect_tolerance_fraction` の 0 の拒否。M4b-1 では `energy_axis` を置かないこと。1次の欄の写像、結果の欄・属性・警告の識別子・CSV の見出し・試験用の口の名前（§4。Codex の表による）。包絡と \(L^1\) の誤差の種（§4.6 の 8）。次数に広げた精度の目標（§7.2）。「値なし」の行の基準（§6.3 の 5）。基準のコードをタグ `v0.3.0` にしたこと。INV の作り方（§6.1）。

**Codex の確認の記録:**

- **1回目（2026-11-27、第1稿に対して）:** 判定は「修正後に凍結」。blocker 4、major 5、minor 1。全文は `work/consult/m4b/prereg/codex_prereg_review.md`。Codex は、表の再現（第1稿の SHA-256 `e255959b…ed80`）、103 行と CI 85 行、現行の M4a での注入の行の1次の値、RB の9条件の1次の成功を確かめた。

| 指摘 | 対応 | 場所 |
| --- | --- | --- |
| 1 [blocker] A3f・AE1・T2e の注入が面積の誤差を2倍にする | 3行とも補間点だけの注入にした。区間の両端は製品の外で作り直してビット単位で照合する | §4.6 の 2、§5.2、表 |
| 2 [blocker] S2 の1次は全次数の誤差の式から外す | 次数 1 は旧欄の転記とし、面積の 11 列の写像を書いた。\(\rho_n\) の式と X22・X24 の対象は 2・3 次 | §4.2、§5.2、§9 |
| 3 [blocker] OI1 の局所の誤差が二分しても消えず、対象のビンがずれている | 親（幅 3 eV 以上）で \(10^{-3}\)、子（1.5 eV 以下）で \(10^{-12}\) の注入にし、対象を `bin[1]` に直した。3次の対照 OI2 を足した | §4.6 の 3、§5.2、§9 |
| 4 [blocker] N19・N20 で面積の (i) が先に超える | `area_epsabs` \(10^{-6}\)。超えた全キーのリストと、代表の選び方を固定し、試験は対象のビンの要素を見る | §4.3、§5.2 |
| 5 [major] 求積点の課金の単位と下限 | 単位を「幅が正の交差区間 \(\times17\)」と固定し、下限を \(527\sum n\) にした。N13 は 4500／4743、N15・C10 は 50,000／10,540 | §3.4、§5.2 |
| 6 [major] D1・D4 の観測の場所 | `on_final_panels`、受け取ったパネルの観測、包絡の種、`_fold_order`、正式の畳み込みの順、3次の対照 OI2 | §3.5、§4.6 |
| 7 [major] 表のキーの分母 | 種類ごとに別の分母で数え、外部の参照を読むキー（第2稿で 9,637、第3稿で 9,668）を分けた。RB1a、A1p、T2c、RB4–6、A1s の台の外の数を直した。行の総数の検査、参照の記録の書式の見本、参照不足の規則 | §5.1、§6.3、§7、表 |
| 8 [major] 境の状態と返す型 | \(N_D=0\) の `"display_support"`（DS0）、`bound_reason` と `tau_per_eV` の3つの場合、下へのあふれ、全軸の残りの全分岐 | §3.1、§3.3、§4.2 |
| 9 [major] 結果の欄と記録の読み方 | 欄の階層、1次の写像、連続区間の最小、解析の入口の \(E_0\) の欄、参照の模型の成分の規約、例外の属性、`invalid_settings` の行（N31–N33）、N18・N21 の使い分け、N27 の順 | §4 |
| 10 [minor] 版・比較幅・プロファイル | 版の欄の規則、精度と被覆の別の判定、CI の行を回しても全件と変異の分母を保つこと | §5.2、§5.3、§7、§8 |
| Q1–Q16、名前の表、変異の経路、運用 | 採った。採らなかった1点（X13 の点の経路の行内の対照）と、変えて採った点は `open_questions.md` の「Codex の確認への返事」 | — |
- **2回目（2026-11-27、第2稿に対して）:** 判定は「修正後に凍結」。blocker 2、major 5、minor 1。1回目の指摘は、1・3・4・5・10 が解消、2・6・7・8・9 が一部と判定された。全文は `work/consult/m4b/prereg/codex_prereg_review2.md`。Codex は、表の再現（第2稿の SHA-256 `a12838d4…7d60`、ハッシュの種 11 と 271）、109／90／19 行、注入の行（A3f、AE1、T2e、N19、N20）の1次の値、DS0・N16 の成立、設計 `69f580f`（§10 の INV と §14 の下限を更新済み）との整合を確かめた。包絡の種（§4.6 の 8）と作業量の下限は採用された。

| 指摘 | 対応 | 場所 |
| --- | --- | --- |
| 1 [blocker] T2c・T4u の \(\tau>0\) を求めるビン | T2c は `bin[0..3]`、T4u は `bin[0..4]`（等号の `bin[4]` を残す）に限り、残りは完全なビンとした | §5.2、表 |
| 2 [blocker] INV の届かない側 | 境より手前は正常の対照（一致）、届く側（またぐビンを含む）で 1 ビン以上の超過を登録。INV1 は 5 eV、INV2 は 11 eV | §5.2、§6.1、表 |
| 3 [major] S2 の次数 1 の `total` と区間の丸め | 面積の内訳に `convention` を足し、再帰の種を別の欄で供給。比の区間と DS0 の \(w\) の区間の作り方（四隅、最後に外向き）を書いた | §4.2、§7.2 |
| 4 [major] 台帳の展開と検査 | DS1 の面積 33 キー。X21 の受け渡しの観測を M4A-o3 に移した。o3・o3s の文を「1回以上」にそろえた。生成の3つの検査（型つきの署名、成果物の供給、変異と具体のキー）を、生成のときの展開で入れた | §5.1、§5.2、§9、生成のスクリプト |
| 5 [major] X13 の点の経路 | A3f に行内の対照を足した（`_fold_order` に渡る節の内訳、絶対の幅 \(10^{-12}\)） | §4.6 の 6、§5.2、§9 |
| 6 [major] A3f の3次の局所の誤差 | 3次の3つの内訳を作り直しにし、AE1 の作り直しと予算の不等式を分けた。分母を出し直した | §5.1、§5.2、§7.2 |
| 7 [major] 派生の量 | 参照の入口を \(Q_h=1\) に限った。高次の感度の誤差の作り方を書いた | §4.2、§4.5 |
| 8 [minor] 総和の順 | 旧版の集計でビット照合すること、「等しい」は実数の式であること、5項は非負の推定であることを書いた | §4.6 の 8 |

- **3回目（2026-11-27、第3稿に対して）:** 判定は「修正後に凍結」。major 3（入力と許容幅を変えずに直せる）。2回目の blocker 2 件は解消、2回目の 4・7 と1回目の 7・9 が一部と判定された。全文は `work/consult/m4b/prereg/codex_prereg_review3.md`。Codex は、表の再現（第3稿の SHA-256 `655e6324…9e48`、ハッシュの種 11 と 271）、台帳の 15,693 件と外部の参照 9,668 件、主入口 109 件の一意性、再利用 12 本、変異の結び先 65 件、A3f の種（0.0368421…）と点の経路の値（0.000603）を確かめた。

| 指摘 | 対応 | 場所 |
| --- | --- | --- |
| 1 [major] X16 に低い KE 側の合計を結ぶ | A1 の台帳に既にある `spectrum.outside.low` を X16 の検出先に足した（行と台帳の件数は増えない） | §9、表 |
| 2 [major] 署名の既定値と、台帳が指せる次数・領域 | 署名は §4.1 の既定値をすべて展開する。台帳と変異の全キーを行（または行内の対照）の次数とビンの範囲に照合する。3つの反例を、`TableError` になる自己検査にした。通し番号だった外部の参照の 52 キーを意味のある名前にした | §5.1、生成のスクリプト |
| 3 [major] 高次の感度に渡す \(\hat r\) の固定の誤差 | 再帰と同じ種 `fixed_multiplier_area_error["term_13a"]` に一本化した。§7.1 の面積の内訳に `convention` を足した | §4.2、§7.1 |

### 凍結の後に実装で確かめる目標

次は、凍結の後に実装で確かめて記録する目標である（Codex の2回目の確認による）。合否の条件を足すものではない。外れたら、値と原因をそのまま記録する。§5.3 の既定の行の移し方で足りない場合や、期待・幅・行・変異を変える場合は、§11 に従う（利用者の了解、旧い表のバイト列、対象のキー、新旧の条件と根拠。比較の後の変更は別の試験として扱う）。

- A1・A1s・A1p・T2c・DS1 の登録したビンと面積で、相対 \(10^{-12}\)、絶対の床はビン \(10^{-13}S_{\rm abs}\)、面積 \(10^{-15}\)。API の終了状態、§7.2 の精度、§7.3 の被覆・過小推定・未確定を別々に記録する。
- T2c・T2e の有限で非ゼロの上界を、同じ \(\tau\) の格子の独立な最小と比べ、下側 \(1-10^{-9}\) と上側 \(1+10^{-6}\) を別々に判定する。定数のパネルの指数の積分、誤差の項の左端の重み、完全なビンを区別する。真の累積の落とした量の被覆と、登録したビンの 10 倍の条件も別に判定する。
- A3f・AE1 の指定した Young の項が参照の値を下回らず、上側の相対の差が \(10^{-9}\) と、明示した入力・丸めの分の内にあること。前提の網が成り立たない場合、参照不足、未実行は合格に入れない。
- CI の増分は、製品、同じプロセスの M4a の対照、基準の子の起動を含む経過 60 秒。超えたら M4A-o3s、C02、RB1 の順に回して毎回測り直す。109 行での CI／全件だけは 90／19、89／20、88／21、87／22（生成が各段階の件数と再利用の向きを検査する）。全件 30 分、各行 120 秒、R-A 5 分、R-B 1コア 30 分、変異 29 件 60 分の枠は保つ。
- X13 の6つの枝、X16 の合計のキー、OI2 の伝播の分と網を作った回数を、実装で観測する。凍結の前に確かめたのは1次の計算と検査の成立の条件までで、高次の精度と 29 の変異の検出の実績は凍結の後の結果である。
- 凍結の時に、本文・表・生成のスクリプトのハッシュとコミットを記録し、独立な参照を固定した後に照合する。

**統括担当の凍結の記録:**（凍結の時に記入する。表の SHA-256、行の数、凍結のコミット）

## 11. 凍結の後の変更

（凍結の時点では空。）凍結の後に期待する状態、許容幅、行、変異を変えるときは、利用者の了解を得て、次を記録する: 変える行とキー、凍結した期待と変えた期待、理由、何が起きたか、原因、確かめたこと。凍結した表はバイト単位で変えずに残し、試験は表を読んだ後に変えた値だけを置き換える（M4a の `POST_FREEZE_EXPECT` と同じ）。参照値を見た後の変更は、旧版を残して別の試験として記録する。

## 付録　凍結の記録（後で記入する）

- 期待する状態の表 `benchmarks/numerics/m4b1/expected_states_m4b1.json` の SHA-256（UTF-8、LF、末尾 LF を含むファイル全体）: `d8d0f892a49c3ed9933b08d6a7c92e5b3aaf05ee5d6c7575050fc164991035e2`
- 生成のスクリプト `expected_states_m4b1.py` の SHA-256: `7c06f0bab21b5f195696fcfe6fbd4bdb0232a7bcf7d9b5f994bd3e2bfbbc248b`
- 凍結の時の設計: `docs/M4B_DESIGN.md` 第5版（コミット `69f580f`）の SHA-256: `6670df1f32db2c77812b7a0e241da226cce6c23753e23929883712374e6dc1ea`
- この文書の SHA-256（自己参照の欄を正規化した値）: `6abbff7acff60bd8fe7bb2bb8bf35d81ec2d199a08350ed22d1d21e337dcef5f`。この行の 64 桁だけを ASCII の `0` 64 個に置き換えたファイルのバイト列（UTF-8、LF、末尾 LF あり）に対して計算する。無加工のファイル全体の SHA-256 は、凍結のコミットとその記録に残す。
- R-A の参照 `references_ra_m4b1.json` の SHA-256、生成のスクリプトと SHA-256、環境、所要時間、キーの数（行ごと）、\(u_r\)/幅 の最大、既知の誤答の記録、「値なし」のキーと当てはめた条件:（参照の担当が、実装の値を見る前に記入。2026-11-27、R-A の担当。実装のコードと結果は見ていない）
  - 参照 `benchmarks/numerics/m4b1/references_ra_m4b1.json`（ASCII、LF、1,836,373 バイト）の SHA-256: `036701007a3d3634b0561a911c1730d54c3ffa81c0d4e3b03f03f09c96ae107a`。
  - 生成のスクリプト `benchmarks/numerics/m4b1/generate_references_ra_m4b1.py` の SHA-256: `4e8beec218c7f3500129738e69b04d41b3587af44b733c512ea19787fc61a0fe`。`xpsloss` を読み込まない。`--check` で再生成とバイト単位で一致する。台帳の具体のキーは、凍結した表の生成のスクリプト（上の SHA-256 を照合してから実行）の展開をそのまま使う。
  - 環境: Windows 11、CPython 3.12.10、mpmath 1.4.1（Python の整数の演算）。厳密な値は `fractions.Fraction`（二進の入力を有理数として読む）、多倍長は 45 桁と 65 桁の2回。値は十進 40 桁。
  - 所要時間: 1プロセスの経過で 約 30 秒（29.5–31.0 秒。上限 5 分）。
  - キーの数: 34 行、3,572 キー（`direct_reference` 3,521、`bound_reference` 38、`error_part` 6、`interval_contains` 7）。表の外部の参照 9,668 から R-B の 6,096 を引いた数に一致。行ごと: A1 264、A1s 114、A1p 114、A2 232、A3a–A3e 各 249、A3f 7、A3z 190、AE1 1、OI1 54、OI2 81、A4a 2、A4c 2、A4b 114、A5 177、I-areas 88、INV1 17、INV2 17、REF1 153、REF2 153、REF3 187、REF4a 27、REF4b 7、REF4c 7、T2a 10、T2c 151、T2e 4、T4v 1、DS1 87、CT1 33、CT7 33。ほかに、台帳の診断のキー3つ（INV1、INV2、T2e）の参照の側の記録を行の `notes` に置いた。
  - \(u_r\)/幅 の最大: \(6.94\times10^{-28}\)（A4a の `order_panels[2].value_at_7eV.term_13a`。幅が相対 \(10^{-13}\) の点の値）。標準の幅の行では \(5\times10^{-34}\) 以下、丸めの水準の行で \(4.7\times10^{-28}\) 以下。多倍長の値（A2、CT7、A5、上界）の \(u_r\) は、2つの作業精度（A5 は求積の次数とパネルの幅も変えた3通り）の差の4倍に、十進への丸めを足した。
  - 「値なし」: 0 キー。当てはめた条件は §6.3 の 1（参照が 0 のキーは \(u_r=0\)。厳密な有理数の 0 だけ）と 2（\(u_r\le\) 幅\(/4\)）。参照の経路の失敗 0、参照不足の行 0。
  - 生成器の自己検査（すべて成立）: 箱形の核で「切断べきの和（部分集合の和の式）＝区分多項式の畳み込み」（A1s、A3a の 2・3 次）、結合則 \((K_1*R)*R=K_1*(R*R)\)（A1s、A3a、A4b。有理数で一致）、全領域の和 \(=Sw_n\)（A2、A4b、A5）、A4a・A4c の \(s=7\) の値が \(\kappa^2\ell\cdot18/19\)、\(\kappa^2\ell\cdot34/33\) に厳密に一致、A5 の \(K*K\) の閉形式が直接の求積と相対 \(10^{-25}\) で一致、INV の逆変換の解（累積を 2 eV の遅れで区分的に解く）が別の関数の次数の和と厳密に一致、REF3 の2次と B1 の2次の差が \(\tfrac12\bar g^{*2}\) に厳密に一致。
  - INV の対照: INV1（境 5 eV）は `bin[9]`–`bin[14]` と `outside.high` で末項を外しても厳密に同じ（正常の対照）。届く側の `bin[0]`–`bin[8]` で末項を外した値が幅の外（例: `bin[0]` で差 \(1.66\times10^{-5}\)、幅 \(3.7\times10^{-9}\)）。INV2（境 11 eV）も届く側に幅の外のビンがある（記録は参照の `notes`）。
  - T2e: `bin[0..3]` のすべてで、誤差の項のない値は誤差の項を入れた参照の 0.826–0.864 倍で、幅（下側 \(1-10^{-9}\)）の外。
  - 既知の誤答（46 件、すべて幅の外）: X01、X02、X04、X05、X06、X07、X08、X09、X10、X12、X14（2通りの読み）、X15、X16、X18（3つの枝×2ビン）、X19、X20、X23、X33（比 0.99785）に加えて、X13（面積の経路と点の経路の各3枝）と X29（T2e の4ビン）。X14 は登録した候補 `bin[5]` が幅の外（粗視化を両方の因子に当てた読みと、1次の側だけに当てた読みの両方）。誤答のキーは、X29 の `bin[1]`–`bin[3]`（登録は `bin[0]`。同じ行の追加の記録）を除いて、すべて表の変異の検出先にある。
- R-B の参照（R-B の担当が、実装の値を見る前に記入。2026-11-27）:
  - `references_rb_m4b1.json` の SHA-256: `ba1d4c8beb7d2d91d713c6355f58b07cced63d261eda5f1c1cf46bc40e01c755`（ASCII、LF、末尾 LF を含むファイル全体。5,033,224 バイト）
  - 生成のスクリプト `benchmarks/numerics/m4b1/generate_references_rb_m4b1.py` の SHA-256: `7ad681e2ed6c3398c570fdb4ddc0958c88ab171ea14835ef315700d48773d6a8`。`--assemble --check` の再生成はバイト単位で同じ。
  - 環境: Python 3.12.10、NumPy 2.5.5、SciPy 1.18.1、Windows-11-10.0.26300-SP0（AMD64）、`xpsloss` 0.4.0.dev0（凍結のコミット `db67d12` の `src`。使ったのは公開の `effective_cross_section_depth_averaged` と `Material`、`Oscillator` だけ）。
  - 方式: 格子 \(h=2^{-k}\) eV（\(k=4..8\)）。値は \(k=5..8\) の3回の Richardson 外挿。点の値は 5 つの組（Al の `k_total`・`k_parallel`・\(Q_h=0\)、Si の `k_total`・`k_parallel`）で 12,545 点または 11,521 点。窓の両端を毎回の呼び出しに入れ、値と誤差が呼び出しによらずビット単位で同じことを確かめた。
  - 所要時間（1コアの経過時間の和）: 点の値 486 秒（Al 71、71、64 秒、Si 139、140 秒）、畳み込みと組み立て 16 秒、計 約 502 秒（8.4 分。上限 30 分の内）。
  - キーの数: 6,096（台帳と一致。RB1–RB8 は各 755 = ビンと両側の外 732、面積 22、\(\hat r\) の区間 1。RB1a は 56 = 面積 22、感度 33、区間 1）。値あり 6,096（うち、その次数の損失が届かない領域の厳密な 0 が 816）。
  - \(f_r\) の最大: \(9.94\times10^{-10}\)（RB1 の `order[3].bin[73].hole`。上限 \(1.25\times10^{-5}\)）。したがって全キーで \(d_{{\rm rel},n}=n\varepsilon+(n-1)\varepsilon_o\)。
  - \(u_r\) の内訳の最大（\(|r|\) に対する比）: (1) 最も細かい2つの \(h\) の外挿の差 \(3.2\times10^{-14}\)、(2) 外挿の回数を変えた差 \(1.1\times10^{-13}\)、(3) 上流の点の誤差 \(1.1\times10^{-10}\)、(4) ビンの一次線の積分の差 \(1.3\times10^{-12}\)、(5) 丸め \(9.5\times10^{-10}\)（\((8M+64)2^{-53}\) と絶対値の和による上側の見積り）。
  - 識別のキー（2次の `total` の最大のビンと両隣）: RB1・RB3・RB7 は `bin[85..87]`、RB2・RB8 は `bin[86..88]`、RB4・RB5 は `bin[84..86]`、RB6 は `bin[83..85]`。すべて値がある。RB1a の必須のキー（2・3 次の `total` の面積と `total` の感度）も値がある。
  - 「値なし」（`no_reference`）: 0 件。参照の経路の失敗: 0 件。参照不足の行: なし。
  - 比較の前の記録（参照のファイルの `pre_comparison_records`）: (a) 2つの運動量の指定の参照の差は、2次の `total` の識別のビンで、2つの比較幅と2つの \(u_r\) の和の Al（RB1 と RB2）で 2,100–7,700 倍、Si（RB4 と RB5）でも和を超える（Al は値のある 99 ビン中 97、Si は 91 ビン中 91 で超える）。(b) X01 の誤った値（2回目以降の核を `total` にする）は、RB1 の2次の `total` の値のある 99 ビンすべて（`bin[19..117]`）で幅を超え、識別のビン `bin[85..87]` では幅の \(3.2\times10^5\)–\(3.8\times10^5\) 倍。3次は 118 領域のうち 117 ビン（`bin[0..116]`）で超える。
  - 自己検査（値には使わない）: 窓の面積 11 列は M4a の R-B の参照（別の経路: 候補の根で分けたガウス則）と相対 \(1.1\times10^{-15}\) 以内で合う。全領域の和と \(Sw_n\)（窓の面積の積）の差は相対 \(8\times10^{-16}\) 以内。
- 完了した比較の表は M4b-1 では作らない（`comparison` は常に `not_matched`）。
