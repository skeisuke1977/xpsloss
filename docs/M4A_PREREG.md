# M4a の事前登録　検証の条件・許容幅・期待する結果

作成日: 2026-11-26。状態: **固定（2026-11-26、実装と参照値の計算の前）。**この文書は、実装と参照値の計算の前に、検証の条件、許容幅、期待する結果を固定する契約である。Codex（実装担当）が第1稿を確かめた（判定は「修正後に凍結」。blocker 2、major 5、minor 2）。第2稿はその指摘を反映した。Codex は第2稿も確かめた（判定は「修正後に凍結」。major 5、minor 2）。第3稿はその修正を反映した。Codex の3回目の確認（判定は「修正後に凍結」。major 1、minor 2）の3点を第4稿に入れた（§10）。4回目の確認はせず、統括担当が確かめて凍結した。

**根拠:**

- [M4a の設計](M4A_DESIGN.md)（第6版）の §3–§6。利用者が §9 の U1–U19 を推奨どおり決定した（2026-11-26）。
- 試走 T1–T16 の報告 `work/consult/m4a/trial/trial_report.md`。**試走の数値は期待値に使わない。**この文書が試走から取ったのは、期待する終了状態、具体的な入力、設計が凍結すると書いた必要なパネル数、時間の見積もりだけである。
- Codex の確認 `work/consult/m4a/prereg/codex_prereg_review.md`（1回目）と `codex_prereg_review2.md`（2回目）、`codex_prereg_review3.md`（3回目）。設計が値や名前を書いていない所は、この確認の案で決めた。

**凍結の後の規則:**

- この文書と期待する状態の表を含むコミットの後に、参照値の計算と実装を始める。
- 比較の後に、許容幅と期待する結果を変えない。変える場合は、この文書の旧版を残し、別の試験として記録する。
- 必須の検査が時間の上限（§8）に収まらなければ、旧版を残して再登録する。
- 不合格なら、層B（セルと一次線の畳み込み）、層A（パネル、解像の尺度、誤差の5項）、駆動と予算の判定、M3a の呼び出し、参照の経路の順に切り分けて記録する。

**範囲:** 原著の式(19)（Simonsen 1997 p.1615）の規約で、表面法線からの出射角 0°、外部一定長 \(L\)、核を \(E_0\) に固定した単一損失の衛星と、有限の損失窓の面積を検証する。式(19)の積分は発散するので（同 p.1618）、窓は必須である。QUEELS-XPS の完全モデルは対象でない。

記号は設計の「表記」に従う。\(K_j\) は M3a の出力（Å⁻¹ eV⁻¹）、\(j\) は 11 列で、順序は `term_13a`、`term_13b`、`term_14`、`term_15_T1`、`term_15_T2`、`term_15_T3`、`term_16_e`、`term_16_h`、`no_hole`、`hole`、`total` である。`no_hole` は 13a、13b、15-T1、15-T3、16-e の和、`hole` は 14、15-T2、16-h の和、`total` は `no_hole`\(+Q_h\cdot\)`hole`。\(\epsilon=2^{-52}\) は binary64 の機械イプシロン、\(I\) は強度の単位である。\(u_r\) は参照値の不確かさ（その量の単位を持つ）で、丸めの単位ではない。

## 1. 入力

数値はすべて十進で与え、API には binary64 に変換して渡す。参照の計算にも同じ binary64 の値を渡す。機械で読める同じ内容が、`benchmarks/numerics/m4a/expected_states_m4a.json` の `materials`、`primaries`、`output_edges` にある（十進と binary64 の16進）。分散は無次元、強度は eV²、他は eV。

### 1.1 材料

| 記号 | 内容 \((E,A,\Gamma,\alpha)\) | 来歴 |
| --- | --- | --- |
| Al | `benchmarks/inputs/al_elfdb_v1.0.bulk.json`（SHA-256 `978931c1651881e6c7dfa739f4d626e4496225bf24bceb815466ff651934eae2`）の材料。\(g=0\)、\((15.6,243.3598,1.5,1.0)\) | ELF データ集 v1.0（CC BY 4.0） |
| Si | `benchmarks/inputs/si_elfdb_v1.0.bulk.json`（`471e513416369eec3dd12302a8e00823113935f78988293205b34f469f7200be`）の材料。\(g=1.12\)、\((10,6.08458,5,0.5)\)、\((14,30.4229,5,0.5)\)、\((16.8,212.9603,3.8,0.5)\) | ELF データ集 v1.0（CC BY 4.0） |
| Al0 | Al の分散だけを 0 にしたもの（CE4） | synthetic |
| SiZ | Si の強度をすべて 0 にしたもの（行 Z1。\(g=1.12\)） | synthetic |
| `synthetic_wide` | \(g=0\)、\((20,200,8,0)\)。CI 用の安い材料 | synthetic |
| N1 | \(g=2\)、\((3,6,0.02,0.5)\) | synthetic（M2-2、M3a と同じ） |
| N2 | \(g=2\)、\((3,6,0.2,0.5)\)。必要なパネル数の表にだけ現れ、行は持たない | synthetic |
| CE1 | \(g=0\)、\((3,3450,0.1,96)\) | synthetic（レビューの反例） |
| CE2 | \(g=0\)、\((1,0.5,0.1,3037.62221442588)\)。分散の16進は `0x1.7bb3e92e3a482p+11` | synthetic（同上） |
| CE3 | \(g=0\)、\((5,0.5,0.1,1.21604938271605)\)。分散の16進は `0x1.374f0329161fcp+0` | synthetic（同上） |
| CE5 | \(g=0\)、\((10,30,1,0)\) と \((10.05,30,1,0)\) | synthetic |
| CE6 | \(g=0\)、\((10,30,1,0)\) と \((10.0000001,30,1,0)\)。第2線の中心の16進は `0x1.40000035afe53p+3` | synthetic |
| P1 | \(g=2\)、\((3,8.1,0.5,0.5)\)（V07 の失敗の対照） | synthetic（M2-2 と同じ） |

### 1.2 共通の条件、一次線、出力の端

**共通の条件:** \(E_0=1000\) eV、\(L=20\) Å（解析の層の行は 2 Å、RA12 だけ 100 Å）、\(Q_h=1\)、`k_total`、出射角 0°。行ごとに変えた点は表に書く。

| 記号 | 一次線（位置・端は eV、強度は \(I\)） |
| --- | --- |
| D1000、D1000x2、D1000zero | δ 1000、強度 1、2、0 |
| D998、D1008、D1200 | δ 998、1008、1200、強度 1 |
| D998_1000 | δ 998 と 1000、強度 1 と 2 |
| D3big | δ 980、1000、1030、強度 1、\(2^{60}\)、2 |
| B_RB3 | ビン端 998、999、1001、積分強度 0.3、0.7 |
| B_RA3 | ビン端 998、999、1002、積分強度 2、1 |

| 記号 | 出力の端（eV） |
| --- | --- |
| E_RB1 | 960、980、984、988、995、999、1000、1001。E_RB1_plus8、E_RB1_plus200 は全点に 8、200 を足したもの |
| E_RB3 | 960、980、984、990、998、1000 |
| E_RA1 | 990、994、995、996.5、998、999、1000、1001 |
| E_RA2、E_RA3、E_RA4 | 996、997、1000／995、996、998、1000／997、998.5、999.5、1000、1001 |
| E_RA5、E_RA6、E_RA6b | 995、996.85、996.87、999／990、1002、1029／990、1002 |
| E_RA7、E_RA7b | 994、995、998、1001／0、0.1、0.2 |
| E_RA9、E_RA10、E_RA12 | 997.98、998／990、1001／997、998 |
| E_Z1、E_L | 990、999、1000、1001／999、999.8、999.95、999.99、1000、1001 |
| E_C04 | `{"min": 960, "max": 1000, "step": 8}`（展開すると 960、968、…、1000 の6点） |

### 1.3 丸めに依る入力（binary64 の16進）

| 使う行 | 値 | 16進 |
| --- | --- | --- |
| N08 | \(a=1.12+5\times10^{-13}\)（上流の除外幅 \(1.12\times10^{-12}\) eV の中） | `0x1.1eb851eb85ab8p+0` |
| N09 | \(a=1.12+2\times10^{-12}\)（除外幅のすぐ外） | `0x1.1eb851eb8751bp+0` |
| N05 | \(h={\rm ulp}(2)=2^{-51}\)、窓 \([2,2+100h]\) | 上端 `0x1.0000000000064p+1` |
| N32b | 窓 \([2,2+256h]\)、\(c=2+64h\)、\(c+h\) | `0x1.0000000000100p+1`、`0x1.0000000000040p+1`、`0x1.0000000000041p+1` |
| N32a | 折れ点 2.5 と `nextafter(2.5, 3)` | `0x1.4000000000000p+1`、`0x1.4000000000001p+1` |
| RA7b | 窓の下端 \({\rm fl}(999.9)\)。\(1000-{\rm fl}(999.9)=0.1+2.27\times10^{-14}>{\rm fl}(0.1)\) | `0x1.f3f3333333333p+9` |
| RA14 | \(k=10^{-10}/0.7753994963\)（binary64 の除算） | `0x1.1b99477b514a2p-33` |
| N16 | 強度 \(2^{901}\)（範囲 \([2^{-900},2^{900}]\) の外） | `0x1.0000000000000p+901` |

### 1.4 既存の関数の比較に使う損失の配列（V01–V06）

GA と GC は作り方が違うので、配列の SHA-256（binary64 のリトルエンディアンのバイト列）を別々に凍結する。表の `v_profiles` に同じ内容がある。

| プロファイル | 配列 | 点数 | SHA-256 |
| --- | --- | ---: | --- |
| CI | GA: `(0.2*np.arange(1, 500))[::4]` | 125 | `dcca1878a047d75b62056854e1bd5344a42ae2c56d570221e0727143e69fa149` |
| CI | GC: 設定 `{"min": 0.2, "max": 99.4, "step": 0.8}` | 125 | `709c4ff510c95d05b14703e2866440aacfc85868a89023a580bff61f90dc4df5` |
| 全件 | GA: `0.2*np.arange(1, 500)` | 499 | `6214a4e35ef20fa9f1c9f5f19a376d6888847f081e5942672eba0ae9b9abf33e` |
| 全件 | GC: 設定 `{"min": 0.2, "max": 99.8, "step": 0.2}` | 499 | `e69297cdd297d6432c4976f071938e792fb93f29f8681cba47da75677a842524` |

CI の GC の上端を 99.4 にしたのは、\((99.8-0.2)/0.8=124.5\) が整数でなく、現在の設定の検査が拒否するためである。

## 2. 既定の精度と上限

**既定（U11。設計 §3.4 の値で凍結する）:**

| 引数 | 記号 | 既定 |
| --- | --- | --- |
| `epsrel` | \(\varepsilon_{\rm rel}\) | \(10^{-6}\) |
| `area_epsabs` | \(\varepsilon_{{\rm abs},A}\) | \(10^{-10}\) Å⁻¹ |
| `spectrum_epsabs_fraction` | \(\eta\)（絶対予算は \(\eta S_{\rm abs}\)） | \(10^{-9}\) |
| `panel_limit` | — | 400 |
| `max_cells` | — | \(2\times10^6\) |
| `kernel_epsabs`、`kernel_epsrel`、`kernel_limit` | M3a の精度 | \(10^{-12}\) Å⁻¹ eV⁻¹、\(10^{-9}\)、200（API の既定） |

- **点の数:** \(n\) 枚の**1つのメッシュ**が持つ異なる損失の点は、共有する端を除いて \(32n+1\)（400 枚で 12,801）。これは累積の評価数の上限ではない。二分の前に評価した親の内部の点は、子の節点にすべては残らず、番兵 \(a\)、\(b\) は段ごとに再び送る。試験は累積の評価数に上限を課さない。段ごとの送信数と累積は `diagnostics.kernel_calls`（§4.4）で読む。
- **予算:** 面積 11 列は \(B=\max(\varepsilon_{{\rm abs},A},\varepsilon_{\rm rel}|A_j|)\)。各ビンと両側の外の \(J_{1,j}\)（3成分）と \(J_0\) は \(B=\max(\eta S_{\rm abs},\varepsilon_{\rm rel}|y|)\)。成功は、判定対象の全量で \(e_y\le B_y\)（推定の基準で、証明された上界ではない。U4、U19 (i)）。

**行ごとに変える設定**（書いていない引数は既定）:

| 行 | 設定 |
| --- | --- |
| CE1–CE6 | `area_epsabs` \(10^{-20}\)（比較の床は \(10^{-18}\) Å⁻¹。§7.2） |
| RA13 | `epsrel` \(10^{-13}\)、`area_epsabs` \(10^{-30}\) |
| RA14 | `panel_limit` 1 |
| N23 | `kernel_limit` 1 |
| N25 | `panel_limit` 9 |
| N26、N27 | `panel_limit` 2、`area_epsabs` \(2\times10^{-4}\)、\(3.5\times10^{-4}\) |
| N28 | `kernel_epsabs` \(10^{-30}\)、`kernel_epsrel` \(10^{-15}\)、`kernel_limit` 200 |
| N30、N36 | `max_cells` 88、89 |
| N31 | `epsrel` \(10^{-17}\) |
| N34 | `area_epsabs` 0.01 |
| N35a、N35b | `panel_limit` 23、24 |

CLI の行（C01–C10c、N29）と V05・V06 は、CONFIG に `"integration": {"epsabs": 1e-12, "epsrel": 1e-9, "limit": 200}` を明示する。省略すると現在の設定の既定 \((10^{-9},10^{-7},200)\) になり、API の既定と違う計算になる。

## 3. 凍結する数え方

期待する状態のいくつかは、次の数え方に依る。数は配列の確保と核の評価の前に決まる。

### 3.1 必要なパネル数（設計 §3.2 (3)(4)）

1. 区切りは、入力の binary64 を厳密な有理数として作る。
2. **刻み:** 低損失は、ギャップなしで \(a\cdot2^k\)、活動のギャップで \(g+(a-g)2^k\)。上端は \(E_0-(E_0-b)2^k\)。比は 2。解析の層は \(E_0\) を持たず、低損失の刻み \(a\cdot2^k\) だけを持つ（RA7b の窓 [\({\rm fl}(999.9)\), 1000.5] も有効である）。
3. **幅の上限:** \(\Gamma_{\min}\)（強度が正の振動子の `width_eV` の最小）。解析の層では、中心のない最小の半幅の 2 倍、または最小の変動の尺度の 2 倍を一律の上限にする。刻みの各区間を \(\lceil\Delta T/w\rceil\)（厳密）に等分する。
4. **許容の規則:** 各パネル \(P\) とすべての候補 \(z_s\) について \({\rm dist}(z_s,P)\ge|P|/2\)。満たさないパネルを binary64 の中点で二分し、満たすまで繰り返す。解析の層で中心と半幅を申告した特異点は、この規則の候補として扱う。
5. 申告した跳び、幅の上限、刻み、許容の規則の二分の境界は併合しない。
6. 必要数が `panel_limit` を超えたら、核を評価する前に `resolution_exceeds_panel_limit`（`required_panels` を付ける）。必要数が上限に等しければ受理する。

**凍結する必要数（初期）:** CE1–CE6 は 14、12、18、74、5、4。N24 は 401（刻みと幅だけの数。候補による二分が起きないことは試走の記録に依る）。RA8a は 41、RA8b は 1,500,002、RA13 は 6、N35 は 24。N25 の `synthetic_wide`・[1,50] は 9（端は 1、2、4、8、16、24、32、38、44、50）。最終の枚数は凍結しない。

**候補の根**（ギャップのない材料。ギャップのある材料は、連成の係数をギャップなしの式で近似する）:

- 下端の根は \(t=1-s/\sqrt{E_0}\)、上端の根は \(t=1+s/\sqrt{E_0}\)（\(s=\sqrt{E_0-z}\)）の多項式で求める。全活動線が無分散なら \(z\) の多項式（\(2n\) 次）で解く。落とすのは、厳密に 0 の最高次の係数だけである。
- 根は元の枝と端に戻して確かめる。族の名前: `bulk`（\(D_i=0\)）、`surface`（\(2-S=0\)）、`epsilon`（\(1-S=0\)）。端は `lower`、`upper`。無分散の固定の極は両端に共通で、`end` を `both` と書く。
- **確認:** 連成の族（`surface`、`epsilon`）は、元の有理式の残差が \(10^{-8}(1+\sum_i|A_i/D_i|)\) 以下。`bulk` は \(|D_i|/\sigma_i\le10^{-8}\)、\(\sigma_i=|W_i|^2+|z|^2+\Gamma_i|z|\)。
- **相殺:** 分母を払った連成の根で、ある \(|D_i|\le10^{-10}\sigma_i\) のものは除く。**重複:** \(|\Delta z|\le10^{-9}(1+|z|)\)。**採用:** \({\rm dist}(z,[a,b])<b-a\)。\(|z|<a/2\) は「原点近傍」の印を付けて残す。
- **根の対応:** 多項式の根と Newton 法で磨いた根を、族・端ごとに多重度を数えて1対1に対応させる。別の根に移った磨きは捨てるが、元の多項式の根は候補に残す。確かめられない根は「未確認」の印（`confirmed: false`）を付けて許容の規則に使い、区分を `computed_candidates_only` に落として警告を出す。

**`resolution_basis`:** `fixed_poles_enumerated`（全活動線が無分散、ギャップなし、全根を確認）、`endpoint_roots_enumerated`（分散あり、ギャップなし、端の全根を確認）、`computed_candidates_only`（ギャップのある材料、または未確認の根がある実行）。活動線のない実行（ゼロの材料の早期の成功）と解析の層では、この欄は `null` で、理由（`zero_response`、`analytic_layer`）を `resolution_basis_reason` に書く。

### 3.2 節点の分離（設計 §3.2 (4)）

各パネルで、核の評価とキャッシュの照会の前に検査する。両端の節点は計算せず、区間の端の値を代入する（\(t_0=\) 左端、\(t_{16}=\) 右端）。内部の 15 の補間点（Chebyshev–Lobatto）と 16 の検証点（第1種 Chebyshev）を中点と半幅から計算して丸めた後、\(t_0<\tau_0<t_1<\dots<\tau_{15}<t_{16}\) が厳密に成り立つことを求める。満たさなければ `panel_width_unresolved`（初期の区切りでは `error`、細分の途中では `not_converged`）。窓の入力の検査は、窓の全体を1枚のパネルとして同じ検査をする（満たさなければ `invalid_window`）。

登録するのは3つの場合だけである: 2 eV からの幅 100 ulp は通らない（N05、窓の拒否）、256 ulp は通る（N32b の窓）、1 ulp は通らない（N32a、N32b の内部のパネル）。設計の「およそ 208 ulp」は、あらゆる入力に当てはめる閾値として登録しない。

### 3.3 セルの数（設計 §3.3 の 8）

出力の領域（各ビンと両側の外）ごとに、窓の全体を、その領域の \(G_y\) の折れ点と跳び、窓の端、パネルの端で区切り、厳密な幅が正の区間をすべて数える。\(G_y=0\) の区間も数える。初期の数が `max_cells` を超えれば `cell_limit_exceeded`（評価 0 回）。細分の各段の前に、その段の後の数を検査し、超えれば `cell_limit_reached`。

**N30・N36 の数:** `synthetic_wide`・[1,50]・D1000・E_RB1 では、初期のパネルは 9 枚である。領域ごとの初期のセルは、低い側の外から順に 10、11、10、10、11、10、9、9、9 で、計 89。88 なら最初の評価の前に入力エラー、89 なら初期は通り、最初の細分の前に未収束になる。

## 4. 実装が守る約束（試験が使う面）

実装担当は参照値を見ない。この節と §5 の表だけで試験を書けるようにする。非公開の名前は Codex の確認の案で決めた。

### 4.1 公開の API（`xpsloss.xps_single_loss`。名前は `xpsloss` からも公開する）

```python
SINGLE_LOSS_MODEL_ID = "simonsen1997_single_loss_spectrum_external_v0"
DeltaPrimary(energies_eV, intensities, intensity_unit="arb. unit", description=None,
             included_broadening=None, truncation=None, intrinsic_satellites=None)
BinnedPrimary(edges_eV, intensities, ...same provenance fields...)
loss_window_area(material, electron_energy_eV, *, loss_window_eV, weight_length_A,
    epsilon_momentum, theta_deg=0.0, hole_charge=1.0, epsrel=1e-6, area_epsabs=1e-10,
    panel_limit=400, kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200,
) -> LossWindowAreaResult
single_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *,
    loss_window_eV, weight_length_A, epsilon_momentum, theta_deg=0.0, hole_charge=1.0,
    epsrel=1e-6, area_epsabs=1e-10, spectrum_epsabs_fraction=1e-9, panel_limit=400,
    max_cells=2_000_000, kernel_epsabs=FIXED_ABS, kernel_epsrel=FIXED_REL, kernel_limit=200,
) -> SingleLossSpectrumResult
```

- 窓、長さ、運動量は必須のキーワード。既存の関数、型、例外、CSV、`run.json` は変えない（V01–V08）。
- 核は公開の `effective_cross_section_depth_averaged` だけを、§4.7 の `_call_kernel` を通して呼ぶ。段ごとにまとめて呼び、毎回 \(a\) と \(b\) を番兵に含める。

**`LossWindowAreaResult` の欄**（すべて float。辞書のキーは 11 列の名前で、上の順序）:

| 欄 | 内容（単位） |
| --- | --- |
| `electron_energy_eV`、`loss_window_eV`（2要素のタプル）、`weight_length_A`、`epsilon_momentum`、`hole_charge`、`theta_deg` | 入力の写し |
| `window_area`、`window_area_error`、`window_area_budget` | 列の名前から float への辞書（Å⁻¹） |
| `single_loss_weight`、`single_loss_weight_error` | \(w_{1,W}=LA_{\rm total}\)（無次元）。誤差は \(Le_{A,\rm total}+\epsilon|w|\) |
| `inverse_window_area_A`、`inverse_window_area_interval_A`、`inverse_reason` | \(1/A_{\rm total}\) と区間（下端, 上端。添字 0、1 が台帳の `lower`、`upper`）\(=[1/(A+e_A),1/(A-e_A)]\)（Å）。定義できるとき理由は `None`。\(A-e_A\le0\) なら値と区間は `None`、理由は `nonpositive_or_unresolved_area`。表せなければ `inverse_not_representable` |
| `cutoff_sensitivity`、`cutoff_sensitivity_error` | \(H_j=LaK_j(a)\)（11 列の辞書、無次元） |
| `diagnostics`、`warnings` | §4.4、§4.5 |

**`SingleLossSpectrumResult` は上に次を加える**（\(n\) は出力のビンの数。配列は `float64`、形 \((n,)\)、出力の端の昇順）:

| 欄 | 内容（単位） |
| --- | --- |
| `output_edges_eV` | `float64`、形 \((n+1,)\) |
| `primary`、`primary_integral`、`primary_abs_integral` | 渡した一次線（δ は `energies_eV` と `intensities`、ビンは `edges_eV` と `intensities`）、\(S\)、\(S_{\rm abs}\)（\(I\)） |
| `zero_loss`、`single_loss_total`、`single_loss_no_hole`、`single_loss_hole`、`spectrum` と、それぞれの `<名前>_error` | ビンごとの配列（\(I\)）。負を切らない |
| `<名前>_per_eV`、`<名前>_per_eV_error` | 密度（\(I\)/eV） |
| `outside["low"]`、`outside["high"]` | 辞書。上の5量の名前と `<名前>_error` が float（\(I\)）。片側の感度は `outside[side]["cutoff_sensitivity"][s][c]` と `outside[side]["cutoff_sensitivity_error"][s][c]`（`s` は `plus`・`minus`、`c` は `total`・`no_hole`・`hole`） |
| `cutoff_sensitivity_bins[s][c]`、`cutoff_sensitivity_bins_error[s][c]` | \(H^\pm_{j,D}\) の配列（形 \((n,)\)、\(I\)） |

- ビンは左閉右開で、最後のビンだけ右端を含む。\(D_-=(-\infty,u_0)\)、\(D_+=(u_N,\infty)\)。台帳の `bin[i]` は各配列の添字 \(i\) である。
- 一次線の台は `diagnostics["primary_support"]` に置く: `min_eV`、`max_eV`（正の強度の点の最小と最大）、`min_minus_E0_eV`、`max_minus_E0_eV`。正の強度がなければ4つとも `None`（Z2）。I3 では2つの差がともに +8 eV になる。
- 感度の `plus` は \(T\to a+0\)、`minus` は \(T\to a-0\) の極限で、\(H^\pm_{j,D}=LaK_j(a\pm0)G_D(a\pm0)\)。最後のビンの右端の扱いと、感度の左右の極限は別の規則である。
- ゼロの材料（全強度が 0）は、共通の検査の後に専用の経路で成功する。核に由来する量は厳密に 0、\(J=J_0\)、逆数は `None`。活動のギャップの検査は適用しない。候補のリストは空、上流の呼び出しは 0 回（Z1）。
- 全一次線の強度が 0 のときは、\(e=B=0\) のキーを成功とし、そのキーの比は `None`、印は `zero_budget_zero_error`。材料の面積は通常どおり判定する（Z2）。

### 4.2 例外と理由

| 型 | 状態 | 理由（`reason`） |
| --- | --- | --- |
| `SingleLossInputError(ValueError)` | `error` | `invalid_window`、`window_not_above_gap`、`invalid_primary`、`invalid_output_edges`、`invalid_settings`、`invalid_schema`、`cell_limit_exceeded`、`resolution_exceeds_panel_limit`、`panel_width_unresolved` |
| 上流の例外 | `error` | M3a の `ValueError`、`DielectricGapEdgeError`、`FixedDepthGapDomainError`。型と公開の属性を保つ |
| `SingleLossConvergenceError(NumericalConvergenceError)` | `not_converged` | `kernel_not_converged`、`kernel_error_exceeds_budget`、`loss_budget_not_met`、`cell_limit_reached`、`panel_width_unresolved`、`nonfinite_output` |

**属性**（該当しないものは `None`）:

| 属性 | 持つ理由 | 内容 |
| --- | --- | --- |
| `reason` | すべて | 上の文字列 |
| `required_panels` | `resolution_exceeds_panel_limit` | 初期の必要数（int） |
| `panel_eV`、`panel_width_ulp` | `panel_width_unresolved`（両方の型） | 分離が成り立たなかったパネルの（左端, 右端）と、その幅の ulp 数 |
| `dominant_part` | `loss_budget_not_met` | `interpolation`、`verification`、`kernel_point_error`、`edge`、`rounding` |
| `node_kind`、`loss_eV`、`index`、`column`、`stage` | `kernel_not_converged` | `sentinel` か `interior`、失敗した損失、その呼び出しの損失の配列の 0 始まりの添字、上流の列、段（初回は 0） |
| `inner_exception`、`inner_reason`、`inner`、`inner_attributes` | `kernel_not_converged` | 上流の型の名前、上流の `reason`、文言、上流の公開の属性の辞書 |
| `diagnostics` | `SingleLossConvergenceError` のすべて | 止まる直前の §4.4 の `diagnostics` と同じ形のスナップショット。結果が返らない行（RA14、N25、N36）でも、`error_budget.outputs["window_area.total"]["estimate"]` などを読める |

**受理と拒否の規則（設計 §5）:**

- 窓: 有限、\(0<a<b<E_0\)、§3.2 の検査。
- 活動のギャップ: \(a\le g\)（厳密な比較）なら `window_not_above_gap`、評価 0 回。除外幅（\(|T-g|\le\delta_g\)）は上流の規則で、最初の M3a の呼び出しが `DielectricGapEdgeError` を送出する（N08）。
- 一次線: δ かビンの一方。強度は非負で有限、正の強度は \([2^{-900},2^{900}]\)、\(S_{\rm abs}\le2^{900}\)。位置とビン端は非負で有限、ビン端は狭義の増加。長さが合うこと、空でないこと。2種類の一次線の組は、評価の前に拒否する（N11）。型の構築での送出も許す。
- 出力の端: 有限、非負、狭義の増加、2点以上。
- 精度: 有限の正、\(\varepsilon_{\rm rel}\ge64\epsilon\)。
- **停止:** 現在の節点の床では止めない。止まるのは、次の段の二分で `panel_limit` を超えるとき（評価の前に検査する）、セルの上限、入力の設定の不正だけである。上限で止まったとき、最悪の出力で点の誤差だけが予算を超えていれば `kernel_error_exceeds_budget`、ほかは `loss_budget_not_met`。

### 4.3 CLI

```powershell
xpsloss xps-loss-window-area CONFIG --loss-window-eV A B --weight-length-A L [--weight-length-source TEXT] --epsilon-momentum M --out DIR
xpsloss xps-single-loss CONFIG --spectrum SPECTRUM.json --loss-window-eV A B --weight-length-A L [--weight-length-source TEXT] --epsilon-momentum M --out DIR
```

- 任意の引数: `--hole-charge`、`--theta-deg`、`--epsrel`、`--area-epsabs`、`--panel-limit`、`--kernel-epsabs`、`--kernel-epsrel`、`--kernel-limit`。`xps-single-loss` だけ `--spectrum-epsabs-fraction` と `--max-cells` を受ける。`--weight-length-source` は両方のコマンドが受ける。
- **精度の優先順:** M3a の精度は、CLI の明示値、CONFIG の `integration`、設定の既定の順。M4a の精度と上限は、CLI の明示値、API の既定の順。出所を `run.json` に記録し、API の再生には `diagnostics.effective_settings` を使う。
- 設定ファイルの `loss_grid_eV` は検証するが使わない（`run.json` に `config_loss_grid_used: false`）。
- **`spectrum.json`:**

```json
{"schema": "xpsloss-single-loss-spectrum-v0", "axis": "kinetic_energy_eV",
 "primary": {"kind": "delta", "energies_eV": [1000.0], "intensities": [1.0],
             "intensity_unit": "arb. unit", "description": null, "included_broadening": null,
             "truncation": null, "intrinsic_satellites": null},
 "output_edges_eV": [960, 980, 984, 988, 995, 999, 1000, 1001]}
```

  - `kind` が `"delta"` なら `energies_eV`、`"binned"` なら `edges_eV` を持ち、`intensities` は共通。片方専用の余分なキー、未知のキー、非有限の値は `invalid_schema` で拒否する。
  - `output_edges_eV` は配列、または `{"min","max","step"}`。後者は \(0\le\) min \(<\) max、区間数が整数、展開の後に狭義の増加。展開した端を昇順の連続した `<f8` のバイト列にして、その SHA-256 を記録する。
- **`window_area.csv`:** 見出し `column,window_area_per_A,window_area_error_per_A,window_area_budget_per_A,cutoff_sensitivity,cutoff_sensitivity_error`。11 行で、列の順序は上のとおり。
- **`single_loss_spectrum.csv`**（スペクトルのコマンドだけ）: 見出しは `region,edge_low_eV,edge_high_eV`、続いて5量の `<名前>,<名前>_error`（`zero_loss`、`single_loss_total`、`single_loss_no_hole`、`single_loss_hole`、`spectrum` の順）、続いて感度の `cutoff_<s>_<c>,cutoff_<s>_<c>_error`（`plus` の `total`、`no_hole`、`hole`、次に `minus` の同じ順）。行は `low`、`0`…`n-1`、`high`。外の無限の側の端は空欄。単位は、端が eV、ほかは \(I\)。
- 浮動小数は既存の CSV と同じ `%.17g` で書く（binary64 に往復できる）。見出しは ASCII。全 CSV を書いた後に `run.json` を完了印として書く。stdout と stderr の JSON は ASCII。
- **失敗:** 終了コード 2。どの失敗でも出力先を作らない。
  - M4a の例外: stderr の JSON に `status`（`error` か `not_converged`）、`error_type`、`reason`、`message`、§4.2 の属性のすべて。非有限の値は文字列にする（C10b では `loss_window_eV` が `[1.0, "nan"]`）。
  - 上流の例外: 既存の欄（`error_type`、`message` と公開の属性）を保つ。
  - 構文エラー（未知の引数、必須の引数の欠落。C09、C10c）: 既存の規則のとおり `status` が `error`、`error_type` が `ValueError`、`message`。`reason` は持たない。
- 試験は、同じプロセスの中で `xpsloss.cli.main(argv)` を呼ぶ（§5.3）。

### 4.4 `run.json`、`diagnostics`、`capabilities`

**`run.json` のパス**（既存の共通の欄 `schema_version`、`status`、`command`、`created_utc`、`elapsed_seconds`、`input_sha256`、`input_file`、`package_version`、`model_id`、`runtime`、`config`、`units`、`constants`、`files`、`warnings` はそのまま持つ）:

| パス | 内容 |
| --- | --- |
| `model_id`、`dependent_model_ids` | `SINGLE_LOSS_MODEL_ID` と、キー `depth_averaged`、`fixed_depth`、`dielectric` を持つモデル ID の辞書 |
| `runtime` | Python、NumPy、SciPy の版 |
| `spectrum_input` | `file`、`sha256`（元のバイト列）、`output_edges_sha256_float64_le`。面積だけの実行では `null` |
| `spectrum_definition` | `loss_window_eV`、`weight`（`length_A`、`length_kind: "external_constant"`、`length_source`（文字列か `null`））、`epsilon_momentum`、`hole_charge`、`theta_deg`、`kernel_energy_eV`、`axis`、`config_loss_grid_used: false` |
| `effective_settings`、`effective_settings_sources` | 計算に使った値と、それぞれの出所（`cli`、`config`、`default`）。キーは実際に API に渡す引数に合わせる。面積は 6 個（`epsrel`、`area_epsabs`、`panel_limit`、`kernel_epsabs`、`kernel_epsrel`、`kernel_limit`）、スペクトルはそれに `spectrum_epsabs_fraction` と `max_cells` を加えた 8 個 |
| `primary` | 一次線の値、\(S\)、\(S_{\rm abs}\)、`support`（`diagnostics.primary_support` と同じ辞書）。面積だけの実行では `primary` 自体が `null` |
| `output_edges_eV` | 展開した端。面積だけの実行では `null` |
| `results` | §4.1 の欄と同じ名前（配列はリスト、`None` は `null`）。`outside` を含む |
| `diagnostics`、`warnings`、`comparison` | 下、§4.5、§4.6 |

**`diagnostics`**（API、`run.json`、例外のスナップショットで同じ名前）:

| 欄 | 内容 |
| --- | --- |
| `resolution_basis`、`resolution_basis_reason` | §3.1 |
| `candidates` | 候補の根のリスト。各要素は `z_re`、`z_im`、`family`、`end`、`residual`、`confirmed`、`near_origin`、`method`（`enumerated` か `search`）、`distance_to_window_eV` |
| `required_panels`、`final_panels`、`stages` | 初期の必要数、最終の枚数、細分の段数 |
| `gradings`、`tolerance_bisected_panels` | 刻みと、許容の規則で二分したパネル |
| `verification` | `rho_max`（生の残差の最大）と `epsilon_max`（検証の点の不確かさの最大）を別々に |
| `error_budget` | `outputs[key]` と `max_ratio`（下） |
| `balance` | 収支の差と丸めの上界（補正しない。合格の根拠にしない） |
| `kernel_calls` | `calls`（呼び出しの回数）、`points_sent_per_stage`（段ごとの送信数のリスト）、`points_sent_total`（累積）、`distinct_points`（異なる座標の数）、`batch_sizes`（呼び出しの順の送信数。番兵の再送を含む。和は `points_sent_total`、長さは `calls`）。I8 の実際の添字は、`_call_kernel` で捕まえた配列から読む |
| `primary_support` | §4.1（スペクトルだけ） |
| `initial_cells`、`cells` | §3.3 の数（スペクトルだけ） |
| `effective_settings` | 計算に使った精度と上限（面積 6 個、スペクトル 8 個） |
| `comparison` | §4.6 と同じ表示 |
| `LbK_b` | \(LbK_j(b)\)（診断だけ） |

- **`error_budget.outputs[key]`:** `key` は §7.1 の名前（`window_area.total`、`bin[3].single_loss_hole`、`outside.low.zero_loss` など）で、判定対象の全量（面積 11、各領域の \(J_0\) と \(J_1\) の3成分）を持つ。各要素は `interpolation`、`verification`、`kernel_point_error`、`edge`、`rounding`、`estimate`（5項の和）、`budget`、`ratio`、`flag`。`flag` は通常 `None`、\(e=B=0\) のキーだけ `zero_budget_zero_error` で `ratio` は `None`。まだ計算していない出力の要素は作らない。初回の上流の失敗（N28、N29）では `outputs` は空の辞書、`max_ratio` は `None` とし、架空の 0 の推定を置かない。
- **`error_budget.max_ratio`:** `{"key", "value"}`。同率なら、面積の列の順、次に領域の順（低い側の外、ビンの昇順、高い側の外）、次に量の順（`zero_loss`、`single_loss_total`、`single_loss_no_hole`、`single_loss_hole`）で先のものを選ぶ。
- **`capabilities`:** `supported` の既存の要素と順はそのままにし、末尾に `"xps_loss_window_area"`、`"xps_single_loss"` をこの順で足す。既存のブロックは等しいままにする。stdout の JSON の `status` は既存のとおり `experimental`（V08）。

### 4.5 警告

`result.warnings` は文言のタプル。試験は `diagnostics["warnings"]`（`{"source", "code", "message"}` のリスト）の識別子で判定する。

| 識別子（出所 `xps_single_loss`） | 条件 |
| --- | --- |
| `window_only`、`weight_not_probability`、`kernel_fixed_at_E0`、`primary_assumed_without_intrinsic_loss` | 常に（最後の1つはスペクトルだけ） |
| `window_over_E0_gt_0p1` | \(b/E_0>0.1\) |
| `window_lower_below_primary_bin_width` | \(a<\) 一次線の最小のビン幅 |
| `E0_outside_primary_support` | \(E_0\) が正の強度の点の最小と最大の間にない。正の強度がなければ適用しない |
| `primary_support_below_window_upper` | 台の下端 \(<b\)。正の強度がなければ適用しない |
| `resolution_basis_computed_candidates_only` | `resolution_basis` が `computed_candidates_only` |
| `outside_matched_scope` | §4.6 の表示に `matched` でない群がある |
| `hole_charge_not_unity` | \(Q_h\ne1\) |

継承した警告は、出所 `xps_depth_averaged`、識別子 `inherited`、上流の文言と、`occurrences`（`[{"stage": i, "loss_eV": [...]}]`）を持つ。同じ出所と文言の警告は1つにまとめ、段と節点の対応は `occurrences` に残す。M3a の共通の警告はそのバッチの全節点に、ギャップ端の個別の警告は該当する節点に対応させる。“Only point values are returned.” は M3a の点値の関数についての文である。

### 4.6 照合の表示

- 表は2つ。(a) 登録表（§5 の表）と、(b) **完了した比較の表**。照合の判定は (b) だけから作る。試走や試験用の表から、本番の (b) を作らない。
- **(b) の場所:** 配布物の中の `xpsloss/data/m4a_completed_comparisons.json` を `importlib.resources` で読む。読む関数は `xpsloss.xps_single_loss._load_comparison_table()` で、試験だけがこれを差し替える。戻り値は辞書で、`version`、`registration_sha256`、`reference_sha256`、`records` を持つ。`records` の各要素は `entry`（`loss_window_area` か `single_loss_spectrum`）、`conditions`（正規化した条件）、`key_status`（`key_status[group][key]` が `{"status", "reason"}`。状態は `pass`・`fail`・`no_reference`・`unresolved`・`not_run`）を持つ。表の `comparison_table_shape` に同じ形がある。
- **`conditions` の形:** キーの名前、入れ子、配列の順は次の例（C02 の正規化した条件。表の `conditions_example_C02` と同じ）のとおりに固定する。数は binary64 の値として比べる。振動子は `[energy_eV, strength_eV2, width_eV, dispersion]` で、設定ファイルの順。面積の入口は `primary` と `output_edges_eV` を持たず、`effective_settings` は面積の 6 個。ビンの一次線は `energies_eV` の代わりに `edges_eV` を持つ。

```json
{"material": {"gap_eV": 0.0, "oscillators": [[20.0, 200.0, 8.0, 0.0]]},
 "electron_energy_eV": 1000.0, "weight_length_A": 20.0, "theta_deg": 0.0,
 "epsilon_momentum": "k_total", "hole_charge": 1.0, "loss_window_eV": [1.0, 50.0],
 "effective_settings": {"epsrel": 1e-06, "area_epsabs": 1e-10, "spectrum_epsabs_fraction": 1e-09,
                        "panel_limit": 400, "max_cells": 2000000, "kernel_epsabs": 1e-12,
                        "kernel_epsrel": 1e-09, "kernel_limit": 200},
 "primary": {"kind": "delta", "energies_eV": [1000.0], "intensities": [1.0]},
 "output_edges_eV": [960.0, 980.0, 984.0, 988.0, 995.0, 999.0, 1000.0, 1001.0]}
```

- **入口ごとの一致:** 正規化した数値列の等号で判定する。条件は、材料の数値、\(E_0\)、\(L\)、角度、運動量の指定、\(Q_h\)、窓、全実効の精度と上限。スペクトルの入口は、それに一次線と出力の端を加える。空白、キーの順、来歴の文字列、長さの出所は比べない（入力の生のハッシュは別に保持する）。
- **表示:** `run.json` の `comparison` と `diagnostics.comparison` に、群 `areas`、`spectrum` ごとの文字列を置く。面積だけの入口は `spectrum` を持たない。詳細は `comparison.key_status[group][key]` に状態と理由を置く。
  1. 同じ入口で全実効条件が一致する完了記録があるとき: 群の必須キーが**すべて** `pass` のときだけ `matched`、それ以外は `not_matched` で確定する。「それ以外」は、キーが記録に欠ける場合と、キーがあって状態が `no_reference`、`fail`、`unresolved`、`not_run` の場合を含む。記録に欠けたキーは、表示の `key_status` で `no_reference` と理由を付ける。
  2. 全条件が一致する記録がなく、核の条件だけが一致して一次線か出力の端が違う記録があるときに限り `kernel_conditions_checked`。
  3. 該当する完了記録がなければ `not_matched`。
- 別の欄に `physical_validity: "not_established"`、`resolution_basis`、計算の `status` を独立に記録する。CE の成功を、一般の面積の入口の `matched` に昇格させない。
- **試験用の表:** `test_table_C02` は C02 の入口と条件を全必須キー `pass` で登録したもの。`test_table_C08_missing_key` は、同じ登録から `spectrum` の必須キー `outside.low.single_loss_total` を欠いたもの。`test_table_C08_no_reference` は、同じキーを状態 `no_reference` で明記したもの（C08 の行内の検査。読み込みの後の同じ判定で、同じ表示になる）。

### 4.7 試験用の口

解析の経路と公開の経路は、同じ細分の駆動と誤差の計算を通す。

1. **`xpsloss.xps_single_loss._call_kernel(material, electron_energy_eV, loss_eV, **kw) -> DepthAveragedResult`。**M3a への唯一の通り道。試験はこれを包んで、呼び出しの回数、損失の配列、実引数（`weight_length_A`、`hole_charge`、`epsilon_momentum`、精度）を記録する。公開の入口の行の「評価 0 回」は、この口が1回も呼ばれないことである。解析の層の行（RA8b、N32a、N32b、N35a など）の「評価 0 回」は、`integrate_window` に渡した `kernel(T, side)` そのものの呼び出しを数える（`_call_kernel` は解析の核を通らないので、この数え方でなければ「評価の前に拒否」を検査できない）。N29 は、ここで初回のバッチから \(a<T<b\) の最初の添字を選び、その損失、添字、列、推定、予算を持つ本物の `DepthAverageConvergenceError` を送出する。正常な C02 が同じ口を1回以上通ること、注入が1回以上発火したことを別に確かめる。
2. **`xpsloss._loss_panels.integrate_window(kernel, loss_window_eV, *, features, weight_length_A, hole_charge=1.0, primary=None, output_edges_eV=None, initial_edges=None, inject=None, epsrel, area_epsabs, spectrum_epsabs_fraction, panel_limit, max_cells) -> (result, panels)`。**
   - `kernel(T, side) -> (values, errors)`。`T` は `float64[n]`、`side` は点ごとの `"left"`・`"right"`・`"interior"` の列、返値は `float64[11, n]`（列の順序は上のとおり）。
   - `features` のキーは `jumps_eV`、`kinks_eV`、`singularities`（`(center_eV, half_width_eV)` の列）、`min_half_width_eV`、`min_variation_scale_eV`。申告のない特徴は契約の外である。
   - `inject(role, T, values, errors) -> (values, errors)`。`role` は `"interpolation"` か `"verification"`。同じ形の `float64[11, n]` を返す。キャッシュの共有で、置き換えた役割を失わない。
   - `initial_edges` は初期の区切りだけを指定し、分離の検査とその後の細分を通す。N26・N27 は `initial_edges=[2, 4, 5]` を渡す（表の `initial_edges_eV`）。
   - `result` は §4.1 と同じ欄と `diagnostics` を持つ（`primary` が `None` なら面積の欄だけ）。
3. **`xpsloss._fold.fold(panels, primary, output_edges_eV, *, multiplier, hole_charge) -> folded`。**層B だけ。核を呼ばない。`folded.outputs[key]` が、領域ごとのキー（§7.1 のスペクトルの群の名前）に対して `value`、誤差の5項（`interpolation`、`verification`、`kernel_point_error`、`edge`、`rounding`）、`estimate` を持つ。予算の設定を受けず、`budget`、`ratio` と予算の判定は持たない。RA1 は 2 の `panels` をそのまま渡し、面積は層A の `result` から読む。RA10 は `folded.outputs["bin[0].single_loss_total"]["value"]`、RA9 は同じ形の要素の `kernel_point_error` を読む。RA9・RA10 は `xpsloss._loss_panels.Panel.from_lobatto(left_eV, right_eV, values, errors)`（各 `float64[11, 17]`、座標の昇順。補間と検証の追加の誤差はこの器では 0）で作ったパネルを渡す。
4. **`xpsloss.xps_single_loss._single_loss_spectrum(..., kernel_batch_size=1)`（I8）。**1 は「1回の呼び出しに、新しい内部の点を1つ」を指し、各呼び出しに重複を除いた番兵 \(a\)、\(b\) を含める。

解析の核の既定の列は `term_13a` で、ほかの項は 0、`no_hole` と `total` は同じ値、`hole` は 0 とする（RA1、RA3、RA10 は表のとおり）。Lobatto の節点の添字は、左端を 0 とする。

## 5. 期待する状態（100 行）

**表:** `benchmarks/numerics/m4a/expected_states_m4a.json`（スキーマ `xpsloss-m4a-expected-states-v4`、ASCII、LF。SHA-256 `250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408`）。生成は `python benchmarks/numerics/m4a/expected_states_m4a.py`（標準ライブラリだけを使い、`xpsloss` を読み込まない）。`--check` で再生成してバイト単位で比べる。

### 5.1 定義と数

- **行**は「新しい入口」（入口、入力、設定の組が他の行と違い、自分の期待する終了状態を持つ実行）を1つ持つ。100 行の（入口、入力、設定）は、正規化して互いに異なる。実行の同一性は、呼び出し、有効な設定、差し替える表で比べ、説明だけの欄（`input_template`、`read_from`）は除く。API の入口は、既定の物理の引数（\(E_0\)、\(L\)、\(Q_h\)、運動量の指定、角度）と数値の設定（面積 6 個、スペクトル 8 個）を展開した実効値で比べる。CLI の入口は、差し替える表と C06 の生の入力の差をキーに残す。
- **行内の検査**（`in_row_checks`）は、行の中で同じ入力・設定の対照の実行（API の再生、乗数だけを変えた層B の再評価、タグの版の実行）と比較を行うもので、行に数えない。
- **再利用**（`reuse`）は、他の行が作った成果物を実際に読むことで、`RB1.result`、`C02.api_replay` のように成果物まで書く。入力のひな型を共有するだけの関係は再利用でない（C04–C08 は各1回の CLI を自分で実行し、C02 の成果物を読まない）。CI の行は、全件だけの行を再利用しない。全件の実行は、表の `execution_order`（再利用の向きに沿った順）で行う。
- **補助の対照**（`auxiliary_controls`）は、行に数えず、合否を行と別に報告する（V07）。
- 行の数の上限は 100。この表は 100 行である。

| 族 | 行 | 入口 | 行数 | CI | 終了状態 | 分母 |
| --- | --- | --- | ---: | ---: | --- | --- |
| RA | RA1–RA7b、RA12 | 解析の層 A・B | 11 | 11 | 成功 | 数値 |
| RA | RA8a、RA13 | 解析の層A | 2 | 2 | 成功 | 数値 |
| RA | RA9、RA10 | 層B だけ | 2 | 2 | 固定パネルで層B が返る | 数値 |
| RA | RA8b、RA14 | 解析の層A | 2 | 2 | 所定の拒否・未収束 | 状態 |
| CE | CE1–CE6 | `loss_window_area` | 6 | 0 | 成功 | 数値（相対）と診断 |
| RB | RB1–RB8 | `single_loss_spectrum` | 8 | 1（RB1） | 成功 | 数値 |
| RB | RB1a–RB4a | `loss_window_area` | 4 | 0 | 成功 | 数値 |
| I | I1a、I1b、I3、I4a、I4b、I8 | `single_loss_spectrum` | 6 | 0 | 成功 | 数値（恒等式と再生） |
| Z | Z1、Z2 | `single_loss_spectrum` | 2 | 1（Z1） | 成功 | 数値（厳密な 0 と恒等式） |
| N | N01–N19、N23–N31、N32a、N32b、N33、N34、N35a、N35b、N36 | 各入口 | 35 | 34（N09 を除く） | 登録した状態（N09、N27、N34、N35b は成功） | 状態 |
| C | C01、C02 | CLI | 2 | 2 | 完了 | 再生と状態 |
| C | C04–C08、C09、C10a–c | CLI | 9 | 9 | 完了（C04–C08）、`error`（C09、C10a–c） | 状態 |
| V | V01–V08 | 既存の API と CLI | 8 | 8 | 登録した状態 | 不変の比較 |
| L | L01–L03 | `single_loss_spectrum` | 3 | 0 | 記録だけ | なし |
| 計 | | | **100** | **72** | | |

- **状態の分母:** 全 100 行（各行に `terminal_state` のキーが1つ）。L の3行は「実行して記録した」ことを状態とし、終了状態は問わない。
- 期待する終了状態の内訳: 成功 49、CLI とハーネスの完了 10、`error` 31、`not_converged` 7、記録だけ 3。
- **数えないもの:** skip、比較 0 件、期待と違う例外は、合格にも検出にも数えない。

**キーの台帳（`check_plan`）:** 各行は、予定するキーを `キー|検査の種類|オラクル` の文字列で、参照値のないいまの段階で列挙する。数値の精度の合格率と、独立な参照の充足率は、この台帳から出す。予定したキーは、結果が `no_reference`、`unresolved`、`not_run` でも台帳に残し、合格の数と区別する。

| 検査の種類 | 内容 | 予定の数 |
| --- | --- | ---: |
| `state` | 行の終了状態（型、`reason`、属性） | 100 |
| `direct_reference` | 独立な参照との精度の比較（§7.2）と推定の三値の検査（§7.3） | 1986 |
| `inverse_interval` | 逆数の区間の2通りの検査（§7.2。§7.3 は重ねない） | 50 |
| `identity` | 恒等式（§7.4） | 132 |
| `exact_zero_or_null` | 厳密な 0、または値・区間・理由の構造 | 169 |
| `replay` | 再生や別の行とのビット単位の一致 | 414 |
| `error_part` | 誤差の内訳の1項と独立な値の比較 | 5 |
| `diagnostic` | 診断の欄（必要数、区分、候補、警告、実引数、照合の表示、CSV の数値の列） | 124 |
| `record` | 記録だけ | 9 |

**必須の数値の比較（`direct_reference` と `inverse_interval`）は 2036 キー**で、内訳は R-A 922、CE 66、R-B 1048 である。

- 面積の群は 25 キー（面積 11、\(w\) 1、逆数の区間の両端 2、\(H\) 11）。スペクトルの群は「ビンの数＋2」の領域ごとに 11 キー（5量と、左右各3成分の感度）。
- R-B: スペクトルの6行（E_RB1 系、7 ビン）が各 124、Si の2行（5 ビン）が各 102、面積の4行が各 25。\(6\times124+2\times102+4\times25=1048\)。
- **CE: 面積 11 列だけで、6 行で 66 キー。**候補、初期の枚数、区分、警告は `diagnostic` で数える。
- R-A: 層 A・B の 11 行が計 869、RA8a と RA13 が各 25、RA10 が 3。RA9 は `error_part` の 3 キーだけ。
- I1a と I4a は独立な参照を持たない。後の行の恒等式に使うだけなので、そのキーを照合済みに数えない。

### 5.2 行の内容

各行の機械で読める定義（`id`、`entry`、`inputs`、`settings`、`expect`、`ci`、`reuse`、`in_row_checks`、`denominator`、`check_plan`）は表にある。ここでは内容と期待を書く。

**R-A（解析の核。\(L=2\) Å、\(Q_h=1\)。核の既定の列は `term_13a`）**

| ID | 核と窓（eV） | 申告 | 一次線、端 | 期待 |
| --- | --- | --- | --- | --- |
| RA1 | [2,5] の矩形。13a 0.04、13b 0.03、15-T1 0.02、15-T3 −0.01、16-e 0.005、14 0.02、15-T2 −0.007、16-h 0.002。窓 [1,8] | 跳び 2、5 | D1000、E_RA1 | 成功。面積は列ごとに \(3c_m\)。衛星の総和 0.60 と \(SLA_{\rm total}=0.60\) を別に比べる。\(J_0\) は [999,1000) に 0、[1000,1001] に 1。行内: 同じ `panels` を `fold` に乗数 3 で渡す（面積は層A の結果で不変、衛星は \(|J_3-1.5J_2|\le4\epsilon|1.5J_2|\)） |
| RA2 | [2,5] で 0.10。窓 [2,5] | 跳び 2、5 | D1000、E_RA2 | 成功。低い側の外 0.2、2つのビンに 0.2 ずつ、\(J_0\) は最後のビン |
| RA3 | 13a \(=0.1\max(0,\min(T-2,(5-T)/2))\)、14 \(=0.02(T-3)\)（[2,5]）。\(Q_h=0.5\)。窓 [1,6] | 跳び 2、5、折れ点 3 | B_RA3、E_RA3 | 成功。面積 0.15（`no_hole`）、0.03（`hole`）、0.165（`total`）。[998,1000] の `hole` は \(-2/225\) |
| RA4a | \(0.03/T\)。窓 [0.02,2] | 窓の中は滑らか | D1000、E_RA4 | 成功。\(A=0.03\ln(b/a)\)、\(H=0.06\) |
| RA4b | 同じ核。窓 [0.2,2] | 同じ | B_RA3、E_RA4 | 成功。行内: \(w({\rm RA4a})-w({\rm RA4b})=0.06\ln10\)（面積の差は \(0.03\ln10\) Å⁻¹） |
| RA5 | \(c/[1+((T-t_0)/\gamma)^2]\)、\(c=0.03\)、\(t_0=3.137\)、\(\gamma=10^{-5}\)。窓 [2,5] | 特異点（3.137、\(10^{-5}\)） | D1000、E_RA5 | 成功。逆正接と対数の閉形式 |
| RA6 | RA2 の核。窓 [2,5] | 跳び 2、5 | D3big、E_RA6 | 成功。低い側 \(J_0=1\)、\(J_1=0.6\)。第1ビン \(2^{60}\)、\(0.6\cdot2^{60}\)。第2ビン \(J_1=1.2\)。高い側 \(J_0=2\) |
| RA6b | 同じ | 同じ | D3big、E_RA6b | 成功。高い側 \(J_1=1.2\)、低い側 0.6 |
| RA7 | [1,5] で 0.10。窓 [2,5] | 跳び 1、5 | D1000、E_RA7 | 成功。\(H^+\) は [995,998) に 0.40、\(H^-\) は [998,1001] に 0.40 |
| RA7b | [999,1001] で 0.10。窓 [\({\rm fl}(999.9)\), 1000.5] | 跳び 999、1001 | D1000、E_RA7b | 成功。\(H^\pm\) は両方 [0.1,0.2) |
| RA8a | RA5 の形で \(\gamma=10^{-6}\)。窓 [2,5]。面積の群だけ | 特異点（3.137、\(10^{-6}\)） | — | 成功。必要数 41。面積は閉形式（およそ \(c\pi\gamma\)） |
| RA8b | 同じ核と窓 | 最小の半幅 \(10^{-6}\) だけ | — | `resolution_exceeds_panel_limit`、`required_panels=1500002`、評価 0 回 |
| RA9 | RA2 の核。`Panel.from_lobatto(2, 5, …)` の1枚。添字 2（左端が 0）の節点だけに、全 11 列の点の誤差 \(e=10^{-6}\)（層B だけ） | — | D1000、E_RA9 | **固定パネルで層B が返る。**駆動の成功の判定は適用しない。`folded.outputs[key]` のビン [997.98, 998] の衛星3成分の `kernel_point_error` が、独立に求めた \(|\omega_2|Le\) と合う（\(\omega_2\) はセル [2, \(1000-{\rm fl}(997.98)\)] での基底の積分で、負。およそ \(-0.00160240\) eV） |
| RA10 | 一定の列: 13a 0.10、14 0.02、`no_hole` 0.10、`hole` 0.02、`total` に識別値 0.50。`Panel.from_lobatto(2, 5, …)`、誤差 0（層B だけ） | — | D1000、E_RA10 | **固定パネルで層B が返る。**`folded.outputs[key]["value"]` で読むビンの `single_loss_total` は 3.0、`no_hole` 0.6、`hole` 0.12（注入した列の積分）。**成分の恒等式（§7.4）と駆動の成功の判定は適用しない** |
| RA12 | \(\kappa[1-\mathcal T_{32}(2T-5)]\)、\(\kappa=4\times10^{-11}\)。窓 [2,3]、\(L=100\) Å | 最小の変動の尺度 0.5（契約を破る故障の注入） | D1000、E_RA12 | 成功。細分が1段以上起こる（面積だけなら初期の1枚で通る）。面積と衛星は閉形式に合う |
| RA13 | \(1/[1+((T-2.25)/0.3)^2]\)。窓 [2,5]。面積の群だけ | 最小の半幅 0.3 | — | 成功。必要数 6。比較は \(d_{\rm rel}=10^{-13}\)（§7.2） |
| RA14 | \(k/[1+((T-2.4)/0.5)^2]\)。窓 [2,3]。Lobatto 点に（値 0、誤差 \(K\)）、検証の点に（値 \(K\)、誤差 \(K\)）を `inject` で入れる | 特異点（2.4、0.5） | — | `loss_budget_not_met`、`dominant_part` は `verification`。行内: 例外の `diagnostics` から読む `window_area.total` の推定が、およそ \(3.5585\times10^{-10}\) Å⁻¹（参照の担当の値と 1 % 以内） |

**CE（レビューの反例。`loss_window_area`、計算は `area_epsabs` \(10^{-20}\)、比較の床は \(10^{-18}\) Å⁻¹）**

| ID | 材料、指定、窓（eV） | 期待（すべて成功） |
| --- | --- | --- |
| CE1 | CE1、`k_total`、[49.5,50.5] | 必要 14、`endpoint_roots_enumerated`。族 `surface`・端 `lower` の候補が、窓までの距離 0.0222 eV（5 % 以内）にある |
| CE2 | CE2、`k_parallel`、[49.5,50.5] | 必要 12、`endpoint_roots_enumerated`。族 `bulk`・端 `lower` の候補が距離 0.02536 eV（1 % 以内）にある |
| CE3 | CE3、`k_total`、[989.5,990.5] | 必要 18。族 `bulk`・端 `lower` の候補が距離 0.00503 eV（1 % 以内）にある。\(b/E_0\) の警告。`resolution_basis` は凍結しないが、`computed_candidates_only` になったときは、その警告が要る |
| CE4 | Al0、`k_total`、[900,999.9] | 必要 74（上端の刻みを外すと 67）。`fixed_poles_enumerated`、\(b/E_0\) の警告 |
| CE5 | CE5、`k_total`、[7,11] | 必要 5。`fixed_poles_enumerated`。\(8.3811-0.5i\) と \(10.0130-0.5i\) のそれぞれから複素の距離 \(10^{-3}\) eV 以内に、別々の `surface` の候補がある |
| CE6 | CE6、`k_total`、[7,11] | 必要 4。\(8.3516-0.5i\) と \(9.9875-0.5i\) のそれぞれから \(10^{-3}\) eV 以内に別々の `surface` の候補があり、第2の根は `confirmed: false`。`bulk` の極（`end: both`）がある。`computed_candidates_only` とその警告 |

CE1–CE4 の数値の一致は、解像の規則の検出にならない。CE の役割は「成功すること」と診断の照合である。規則の検出は RA8b、N32b、N24 の状態と、CE の必要数の照合で行う。

**R-B（実材料）**

| ID | 条件（共通の条件から変えた点） | 一次線、端 | 入口 |
| --- | --- | --- | --- |
| RB1、RB2 | Al、[1,50]、`k_total`・`k_parallel` | D1000、E_RB1 | スペクトル（RB1a、RB2a は面積） |
| RB3、RB4 | Si、[5,50]、`k_total`・`k_parallel` | B_RB3、E_RB3 | スペクトル（RB3a、RB4a は面積） |
| RB5 | Al、[1,50]、\(L=5\) Å | D1000、E_RB1 | スペクトル |
| RB6、RB7 | Al、[1,50]、\(Q_h=0\)、0.5 | D1000、E_RB1 | スペクトル |
| RB8 | Al、[1,50]、\(E_0=1200\) eV | D1200、E_RB1_plus200 | スペクトル |

すべて成功。Si の行は `computed_candidates_only` の警告、RB6・RB7 は \(Q_h\) の警告を持つ。`_call_kernel` に渡った `weight_length_A`、`hole_charge`、`epsilon_momentum` も確かめる。**識別のビン**（M06 の検出先）: Al は `bin[1]`（[980,984)）と `bin[2]`（[984,988)）、Si は `bin[0]`（[960,980)）と `bin[2]`（[984,990)）の `single_loss_total`。全件の実行では、RB1 の照合の表示を行内で確かめる。

**恒等式・ゼロ・低損失**

| ID | 内容 | 期待 |
| --- | --- | --- |
| I1a、I1b | Al、[1,50]、E_RB1。一次線は D998、D998_1000 | 成功。I1a は `E0_outside_primary_support` の警告。行内（I1b）: \(J(F_1+2F_2)=J(F_1)+2J(F_2)\)（`I1a.result` と `RB1.result` を読む。§7.4） |
| I3 | RB1 の一次線と端を +8 eV（D1008、E_RB1_plus8）。\(E_0\) は 1000 eV のまま | 成功。数値と誤差が `RB1.result` とビット単位で同じ。端と一次線の位置の欄、\(E_0\) との差の欄、`E0_outside_primary_support` の警告は変わる |
| I4a、I4b | Al、D1000、E_RB1、窓 [1,17] と [17,50] | 成功。行内（I4b）: `RB1.result` との加法性（§7.4） |
| I8 | `synthetic_wide`、[1,50]、D1000、E_RB1、`kernel_batch_size=1`。実効の設定は C02 と同じ（M3a の精度 \(10^{-12}\)、\(10^{-9}\)、200） | 成功。数値と誤差のバイト列が `C02.api_replay` と同じ。バッチの数と添字は `kernel_calls` で別に確かめる |
| Z1 | SiZ、D1000、E_Z1、窓 [0.5,2]（ギャップをまたぐ） | 成功。核に由来する量は厳密に 0、逆数は値と区間が `None` で理由が `nonpositive_or_unresolved_area`、\(J=J_0\)。`resolution_basis` は `null`（理由 `zero_response`）、候補は空、上流の呼び出しは 0 回 |
| Z2 | Al、[1,50]、D1000zero、E_RB1 | 成功。\(J_0=J_1=J=0\)。正の強度の台と \(E_0\) との差は `null`、台に依る警告は出ない。\(e=B=0\) のキーだけ比が `null` と `zero_budget_zero_error`。面積は通常どおり判定する。行内: 面積が `RB1a.result` と、2つの誤差の和の内で合う |
| L01–L03 | Al、\(b=50\)、\(a=0.2\)、0.05、0.01、D1000、E_L | 記録だけ（合否なし）。終了状態、\(H_j\)、成分、誤差の5項を記録する |

**拒否・未収束**

| ID | 入口と入力 | 期待 |
| --- | --- | --- |
| N01–N05 | 面積、Al。窓 [0,50]、[50,50]、[1,1000]、[nan,50]、\([2,2+100h]\) | `invalid_window`、評価 0 回 |
| N06、N07 | 面積、Si。窓 [0.5,1.0]、[1.0,5] | `window_not_above_gap`、評価 0 回 |
| N08 | 面積、Si。窓 [`0x1.1eb851eb85ab8p+0`, 5] | 最初の `_call_kernel` が `DielectricGapEdgeError`。型と公開の属性を保つ |
| N09 | 面積、Si。窓 [`0x1.1eb851eb8751bp+0`, 5] | 成功（全件だけ） |
| N10–N16 | スペクトル、Al、[1,50]、E_RB1。一次線: 空／δ とビンの組／長さの不一致（位置 2、強度 1）／強度 −1／強度 `inf`／ビン端 998、1001、999／強度 \(2^{901}\) | `invalid_primary`、評価 0 回 |
| N17–N19 | スペクトル、Al、[1,50]、D1000。端: [1000]／[960,980,980]／[−1,0,1] | `invalid_output_edges`、評価 0 回 |
| N23 | 面積、`synthetic_wide`、[1,50]、`kernel_limit` 1 | 上流の `ValueError`（文言に `integration limit too small`） |
| N24 | 面積、N1、[2.1,10]、`k_total` | `resolution_exceeds_panel_limit`、`required_panels=401`、評価 0 回 |
| N25 | 面積、`synthetic_wide`、[1,50]、`panel_limit` 9 | 初期の 9 枚を評価した後、次の段の二分の前の検査で `loss_budget_not_met`（`dominant_part` は `interpolation`） |
| N26、N27 | 解析の層A。RA2 の核、窓 [2,5]、`initial_edges` 2、4、5、`panel_limit` 2。`inject` で補間点だけに点の誤差 \(10^{-4}\)、検証点は誤差 0。`area_epsabs` \(2\times10^{-4}\)、\(3.5\times10^{-4}\) | N26: `kernel_error_exceeds_budget`。N27: 成功。行内（N27）: `window_area.total` の `kernel_point_error` が \((2+1)\times10^{-4}\) Å⁻¹（相対 \(10^{-12}\) 以内） |
| N28 | 面積、`synthetic_wide`、窓 [5,6]、`kernel_epsabs` \(10^{-30}\)、`kernel_epsrel` \(10^{-15}\) | `kernel_not_converged`。`node_kind` は `sentinel`、`loss_eV` 5.0、`index` 0、`stage` 0、`inner_exception` は `DepthAverageConvergenceError`、`inner_reason` は `panel_floor_exceeds_budget`。例外の `diagnostics.error_budget.outputs` は空、`max_ratio` は `None` |
| N29 | CLI（スペクトル）。CONFIG と argv は C02 と同じ（`synthetic_wide`、`loss_grid_eV` `{"min": 1, "max": 100, "step": 1}`、\(E_0=1000\)、`integration` 明示、`--loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total`、D1000、E_RB1）。§4.7 の 1 の注入（`budget_not_met`） | `not_converged`、`kernel_not_converged`。`node_kind` は `interior`、添字・損失・列・段と `inner_attributes` は注入した例外のもの（注入のときの辞書と stderr を比べる）。`outputs` は空、`max_ratio` は `None`。終了コード 2、出力先なし |
| N30 | スペクトル、`synthetic_wide`、[1,50]、D1000、E_RB1、`max_cells` 88 | `cell_limit_exceeded`、評価 0 回 |
| N31 | 面積、Al、[1,50]、`epsrel` \(10^{-17}\) | `invalid_settings`、評価 0 回 |
| N32a | 解析の層A。[2,3] の矩形、折れ点 2.5 と `nextafter(2.5,3)` を申告 | `panel_width_unresolved`（`error`）、評価 0 回 |
| N32b | 解析の層A。窓 \([2,2+256h]\)、\([c,c+h]\) で \(4K_0u(1-u)\)（\(K_0=10^6\)）、折れ点 \(c\)、\(c+h\) と尺度 \(h/2\) を申告 | 窓の検査を通った後、`panel_width_unresolved`（`error`）、`panel_width_ulp` 1、評価 0 回 |
| N33 | 解析の層A。[2,3]、核が `inf` を返す | `nonfinite_output` |
| N34 | 解析の層A。[2,3]、核は 1、全点に点の誤差 \(\exp[-((T-2.5)/0.001)^2]\)、`area_epsabs` 0.01 | 成功。親のパネルで止めず、細分が1段以上起こる |
| N35a、N35b | 解析の層A。RA2 の核、窓 [2,5]、最小の半幅 0.0625（必要数 24）、`panel_limit` 23、24 | `resolution_exceeds_panel_limit`（`required_panels=24`）／成功 |
| N36 | N30 と同じ入力で `max_cells` 89 | `cell_limit_reached`（最初の細分の前の検査）。`initial_cells` 89 |

**CLI（`synthetic_wide`、\(E_0=1000\) eV、窓 [1,50]、\(L=20\) Å、`k_total`。CONFIG は `integration` を明示）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| C01 | 面積。`--weight-length-source "design example"` | 完了、終了 0。`run.json` の `primary` は `null`、`spectrum_definition.config_loss_grid_used` は `false`。行内: `effective_settings`（6 個）からの API の再生と面積の群の 25 キーがビット単位で一致。`results` と、`window_area.csv` の全数値の列（誤差と予算を含む）の両方を API と照合する |
| C02 | スペクトル。D1000、E_RB1（配列）。表は `test_table_C02` | 完了。行内: `effective_settings`（8 個）からの API の再生と 124 キーがビット単位で一致。`results` と、2つの CSV の全数値の列の両方を API と照合する。照合の表示は両群 `matched` |
| C04 | 出力の端を E_C04（`min`・`max`・`step`）にする | 完了。両群 `kernel_conditions_checked` |
| C05 | `--hole-charge 0.5` | 完了。両群 `not_matched` |
| C06 | `spectrum.json` の空白・キーの順・来歴の文字列と、長さの出所を変える | 完了。両群 `matched` |
| C07 | 強度を 2 にする | 完了。両群 `kernel_conditions_checked` |
| C08 | C02 と同じ入力。表は `test_table_C08_missing_key` | 完了。`areas` は `matched`、`spectrum` は `not_matched`、`outside.low.single_loss_total` は `no_reference`。行内: `test_table_C08_no_reference` でも同じ表示になる |
| C09 | 面積のコマンドに `--max-cells 10` | `error`、`error_type` は `ValueError`、`reason` なし、終了 2、出力先なし |
| C10a | `spectrum.json` の `axis` が `"binding_energy_eV"` | `error`、`SingleLossInputError`、`invalid_schema`、終了 2、出力先なし |
| C10b | `--loss-window-eV 1 nan` | `error`、`invalid_window`、終了 2、出力先なし。stderr の `loss_window_eV` は `[1.0, "nan"]` |
| C10c | `--epsilon-momentum` を与えない | `error`、`error_type` は `ValueError`、`reason` なし、終了 2、出力先なし |

C04–C08 は、表示と警告 `outside_matched_scope` の有無の対応も確かめる（両群 `matched` の C06 では出ない）。

**既存の API と CLI の不変（タグ `v0.2.0`、コミット `023820cfbc95547a1cc90c11d1f2bcdf19400806` との比較。`k_total`、\(Q_h=1\)、精度 \((10^{-12},10^{-9},200)\) を明示）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| V01、V02 | `effective_cross_section_fixed_depth`、Al・Si、深さ 10 Å、GA | 成功。全欄がビット単位で同じ |
| V03、V04 | `effective_cross_section_depth_averaged`、Al・Si、20 Å、GA | 同上 |
| V05、V06 | CLI `xps-fixed-depth`（深さ 10 Å）、`xps-depth-averaged`（20 Å）、Al、GC | 完了。CSV、`run.json`、stdout の JSON が同じ（stdout の `output_dir` は下の規則） |
| V07 | `effective_cross_section_depth_averaged`、P1、20 Å、1.999 eV | `FixedDepthGapDomainError`、`surface_pole_not_excluded`。型、文言、公開の属性がタグと同じ |
| V08 | `xpsloss capabilities` | 終了 0（表の `complete` はハーネスの正常終了の意味）。stdout の `status` は `experimental`。`version` は各版の `__version__` と照合し、比較から外す。`supported` の旧い接頭部、既存のブロック、追加の2つの名前を別々に検査する |

- GA と GC は §1.4 のプロファイルに従う（CI は 125 点、全件は 499 点）。各環境で V の 8 行を保ち、実行したプロファイル、実際の点数、再生、skip を記録する。
- **V07 の補助の対照**（行に数えず、別に報告する）: (i) P1 の固定深さ（10 Å、1.999 eV、`FixedDepthGapDomainError`）。既存の `tests/test_m3a_migration.py`（M3a の V14(b) の失敗の対照）に当たる。(ii) Si の深さ平均（20 Å、`np.r_[GA_full[:10], 1.12]`（全件の GA の最初の 10 点と 1.12。M3a の GA10G と同じ）、`DielectricGapEdgeError`、添字 [10]）。既存の `tests/test_m3a_expected_states.py`（M3a の RJ1）に当たる。
- **タグの側の実行:** `v0.2.0^{commit}` の SHA を上の値に固定し、`git archive` で取り出した版を、隔離した1つの子プロセスで V の全行に使う。読み込んだモジュールのパスを確かめる。時刻（`created_utc`）、所要時間（`elapsed_seconds`）、絶対パス（`input_file`）の欄は比較から外す。V05・V06 の stdout の `output_dir` は絶対パスなので、各版で実際の出力先と一致することを確かめてから比較から外す。版番号（`run.json` の `package_version`、`capabilities` の `version`）は各版の `__version__` と別に照合する。CI で履歴がなければ飛ばさずに失敗する。展開した sdist では理由を書いて skip する（M3a の V14 と同じ）。

### 5.3 CI と全件

- **CI（72 行）:** R-A の 17 行、N の 34 行、Z1、C の 11 行、V の 8 行、RB1。**全件だけ（28 行）:** CE1–CE6、RB2–RB8、RB1a–RB4a、I の 6 行、Z2、N09、L01–L03。
- **CI の規則:**
  1. N28 は `synthetic_wide`・窓 [5,6] で行う。
  2. N09（除外幅のすぐ外）は全件だけにする。
  3. CLI の行（C の 11 行と N29）は、同じプロセスの中で `main()` を呼ぶ。
  4. タグ `v0.2.0` の側の実行は、V の全行をまとめた1つの子プロセスで行う。
  5. V は CI で 125 点のプロファイル、全件で 499 点のプロファイルを使う。
- **60 秒を超えたときの順（表の `ci_reduction_order`）:** 代表の各 CI の環境で増分の経過時間を測る。超えたら、まず C05、次に C07 を全件だけに回す（71、70 行）。C02、C04、C06、C08、N29 は残す。それでも超えれば、旧版を残して再登録する。回した行と測った時間を記録する。
- 全件の記録は、統括担当が全 100 行を `execution_order` で1回実行し、結果、所要時間、分母ごとの件数を検証記録に書く。

## 6. 参照値（実装担当とは別の担当）

参照の担当は、実装のコードと結果を見ずに作る。実装担当は参照値を見ない。参照値は凍結して SHA-256 を付録に記録してから、実装と比べる。値は十進 30 桁以上の文字列で保存し、\(u_r\) に十進への丸めの幅を含める。参照は、表の `check_plan` のキー（`direct_reference`、`inverse_interval`、`error_part`）に対して作る。

### 6.1 R-A（解析の参照）

- **パッケージ `xpsloss` を読み込まない。**二進の入力を厳密な有理数として読み、閉形式を mpmath（40 桁以上）か有理数で評価する。実装のパネル、セル、重なりの関数を使わない。
- δ 線 \((p,s)\) のビン \([u,v)\) の衛星は \(Ls\int_{W\cap[p-v,\,p-u]}K_j\,dT\)。ビンの一次線は、同じ式を一次線の区間で積分する。両側の外は、その領域に入る部分だけを直接に積分し、全体との差を使わない。
- \(u_r\) は、作業精度を変えた差と、十進への丸めの和。
- RA9 の \(\omega_2\)（[2,5] の 17 点の Lobatto 基底のうち添字 2 のものを、セルで積分した値）と、RA14 の推定の値は、参照の担当が自分で計算して凍結する。
- **RA13 の凍結前の確認:** 参照の担当が、自分の試作で、正常版が §7.2 の幅を通り、初期の 6 枚で打ち切った版が落ちることを確かめて記録する（M14a）。
- **既知の誤答:** §9 の変異に当たる誤った値（例えば RA1 の総和 0.30、1.0、RA3 の 0.18、RA10 の 0.72）を参照の側で作り、§7.2 の幅の外にあることを記録する。

### 6.2 R-B（実材料の通し経路）と CE の参照

- **方式:** 一次線のエネルギー \(E'\) を外側、損失 \(T\) を内側に取る。内側は \(T\in W\cap[E'-v,E'-u]\)（低い側の外は \(T>E'-u_0\)、高い側の外は \(T<E'-u_N\)）。内側の \(T\) と外側の一次線の区間を、端の切り替わりと共鳴で区切った Gauss–Legendre で積分する。δ はその位置で評価する。面積は窓の全体で別に積分する。
- **独立性:** 公開の `effective_cross_section_depth_averaged` だけを共有する。実装のパネル、セル、重なりの関数、節点を借りない。**候補の根（区切りに使う共鳴）は、実装とも試走の試作とも別に計算する。**`total` は M3a の `total` を直接に積分し、ほかの2列から作り直した値だけを参照にしない。窓の両端を含めて M3a の定義域を検査する。
- **\(u_r\) の内訳**（キーごとに別々に記録して足す）: (1) 内側だけ次数を変えた差、(2) 外側だけ次数を変えた差、(3) 区切りを局所で変えた差、(4) 上流の誤差（M3a の列ごとの点の誤差を、求積の重みの絶対値で伝えたもの）、(5) 丸め。\(r\ne0\) のとき \(f_r=u_r/|r|\)。
- **CE の参照:** 同じ型の経路で、面積 11 列だけを作る（66 キー）。候補のまわりを等比に区切る。
- **比較の前に記録すること:** 2つの運動量の指定の独立な参照の差が、識別のビン（§5.2）で2つの比較幅の和を超えること。M04 の誤った実装の値が、**RB5** の有効な参照のキーで幅を超えること。RB1 は元から 20 Å なので、M04 の後も呼び出しと計算値が正常版と一致する（正常の対照。幅の超過を求めない）。M05a（RB6、RB7）の誤った実装の値が、有効な参照のキーで幅を超えること。

### 6.3 「値なし」の基準（実装の値を見る前に凍結する）

キーごとに、次の順で決める。どの条件で「値なし」にしたかを記録する。

1. **参照が 0 のキー（\(r=0\)。例えば RB1 の `outside.high.zero_loss`）:** \(f_r\) は定義しない。床 \(d_{\rm abs}\) が有限で正なら、\(u_r\le d_{\rm abs}/4\) を確かめて値を持つ。比較幅は \(d_{\rm abs}\) そのものである（§7.2）。床も 0 の厳密な 0（`exact_zero_or_null`）は、\(u_r=0\) と \(y=0\) を求める。
2. **逆数が定義できないとき:** 値、区間、理由の構造を比べる（数値の比較にしない）。
3. **\(r\ne0\) の R-A:** 細かくしても \(u_r\le\)（§7.2 の幅）\(/4\) を満たさないキーは `no_reference`。
4. **\(r\ne0\) の R-B と CE:** \(d_{\rm rel}=\max(8f_r,\varepsilon_{\rm rel})\) なので、幅の 1/4 の条件は自動的に通る。有効な精度を縛るのは \(d_{\rm rel}\le10^{-4}\)（すなわち \(f_r\le1.25\times10^{-5}\)）で、超えるキーは `no_reference`。
5. 参照の経路が失敗したキー（上流の未収束を含む）は「参照の経路の失敗」として別に数える。

`no_reference` のキーは、合格にも不合格にも数えず、台帳には残す。比較の後に「照合未確定」（`unresolved`）になったキーは、比較する担当から参照の担当に、キーと必要な精度だけを渡す。参照の担当は、細かくする規則、旧値、新しい値、新旧の SHA-256 を記録する。未確定を合格の件数に含めない。

**別々に数えるもの:** 合格、不合格、`no_reference`、`unresolved`、`not_run`、実装の予期しない未収束、上流の失敗、参照の経路の失敗。

## 7. 比較のキーと許容幅

### 7.1 キー

- **面積の群（25）:** `window_area.<列>`（11）、`single_loss_weight`、`inverse_window_area_interval_A.lower`・`.upper`、`cutoff_sensitivity.<列>`（11）。
- **スペクトルの群（領域ごとに 11）:** 領域は `outside.low`、`bin[0]`…`bin[n-1]`、`outside.high`。量は `zero_loss`、`single_loss_total`、`single_loss_no_hole`、`single_loss_hole`、`spectrum` と、`cutoff_sensitivity.<plus|minus>.<total|no_hole|hole>`。
- 行ごとの必須キーは表の `check_plan` が決める。合否は「行、キー」の単位で数える。`hole`、`no_hole`、`total` は別に判定する。密度（`*_per_eV`）は、ビンの値と幅からの再計算と \(2\epsilon\) 以内で一致することだけを確かめる。

### 7.2 精度の比較

\[
|y-r|+u_r\le d_{\rm abs}+d_{\rm rel}|r|
\]

比較幅は実装の \(e_y\) と独立に決める。**\(r=0\) のキーでは、どの参照でも比較幅を \(d_{\rm abs}\) とし、\(u_r\le d_{\rm abs}/4\) のとき \(|y|+u_r\le d_{\rm abs}\) と §7.3 を検査する**（\(f_r\) と \(d_{\rm rel}\) を使わない）。

| 参照 | \(d_{\rm rel}\)（\(r\ne0\)） |
| --- | --- |
| R-A | 実行の `epsrel`（既定 \(10^{-6}\)。RA13 は \(10^{-13}\)） |
| R-B、CE | \(\max(8f_r,\varepsilon_{\rm rel})\)。最大は \(10^{-4}\)（§6.3） |

| キー | 単位 | \(d_{\rm abs}\) |
| --- | --- | --- |
| 面積 | Å⁻¹ | 実行の `area_epsabs`（既定 \(10^{-10}\)、RA13 は \(10^{-30}\)）。CE の比較は \(10^{-18}\) |
| \(w\) | 無次元 | \(L\times\) 面積の床 |
| \(H_j\) | 無次元 | \(La\times10^{-12}\) |
| ビンと外の5量 | \(I\) | \(\eta S_{\rm abs}\) |
| \(H^\pm_{j,D}\) | \(I\) | \(La\times10^{-12}\times S_{\rm abs}\) |
| 逆数の区間の両端 | Å | 下の2通りの検査 |

- **置き換えの幅（RA6、RA6b）:** 強度を掛ける床（ビンと外の5量と \(H^\pm\)）に限り、\(S_{\rm abs}\) を「その領域にゼロ損失か衛星が届く一次線の強度の和」に置き換える。和は入力の台から決める。低い側の外は 1、RA6 の第2ビンと高い側の外、RA6b の高い側の外は 2、第1ビン [990,1002) は \(2^{60}\)（届くのは 1000 eV の線だけ。980 eV の線は衛星も低い側の外に入り、1030 eV の線の衛星は [1025,1028] eV に入る）。駆動の予算に使う全強度は \(S_{\rm abs}=2^{60}+3\) で、置き換えの強度と書き分ける（この入力では2つの binary64 は同じ値になる）。既定では、衛星 0.6、1.2 の全体の幅は \(6.01\times10^{-7}\)、\(1.202\times10^{-6}\,I\) になる。面積と \(H_j\) には、その単位の床を使う。
  - これは実装の成功の契約より狭い**独立な精度の比較**である。駆動の予算 \(B=\max(\eta S_{\rm abs},\varepsilon_{\rm rel}|y|)\) は U4・U19 のとおりで、\(S_{\rm abs}=2^{60}+3\) では約 \(1.15\times10^9\,I\) になり、裾を丸ごと失っても予算の判定だけでは見つからない。精度の比較と、実際に報告された誤差による三値の検査（§7.3）を別々に記録する。
- **RA3 のビンの `hole`:** \(3\times10^{-9}+10^{-6}|r_{\rm hole}|\)（\(I\)）。識別のビンではおよそ \(1.19\times10^{-8}\,I\)。
- **逆数の区間（2通り。両方を判定する）:** (i) 参照との整合: \(r\)、\(u_r\) は同じ行の `window_area.total` の参照の値と不確かさである（逆数だけの参照は作らない）。参照の \([1/(r+u_r),1/(r-u_r)]\) が実装の区間に含まれれば合格、交わらなければ不合格、ほかは照合未確定。(ii) 実装が報告した \(A\)、\(e_A\) を厳密な有理数として、\(1/(A+e_A)\) を下へ、\(1/(A-e_A)\) を上へ丸めた binary64 と、実装の区間の両端がビット単位で一致する。`inverse_window_area_A` は \(1/A\) の最近接の丸めと一致する。端ごとに架空の \(e_y\) を作って §7.3 を重ねない。
- **RA13:** \(|y-r|+u_r\le10^{-30}+10^{-13}|r|\)。
- **`error_part` のキー:** RA9 は \(|\omega_2|Le\) と相対 \(10^{-9}\)、RA14 は 1 %、N27 は相対 \(10^{-12}\)。

### 7.3 推定の検査（三値）

`direct_reference` のキーで、\(d=|y-r|\)、\(e=e_y+f_{\rm fp}\)（\(f_{\rm fp}=4\epsilon\max(|y|,|r|)\)）として判定する。\(e_y\) は実装が実際に報告した誤差である。

- \(d+u_r\le e\): **推定が覆う**。
- \(d-u_r>e\): **過小推定**（不合格）。
- その間: **照合未確定**（§6.3 の手順で参照を細かくする）。

精度の比較（§7.2）と推定の検査は別の検査で、両方を記録する。公開の入口と `integrate_window` で成功した行の判定対象の量では、\(e_y\le B_y\) も確かめる。層B だけの行（RA9、RA10）には、\(e_y\le B_y\) を課さない。

### 7.4 行内の恒等式の幅

| 検査 | 幅 |
| --- | --- |
| I1b の線形性 | 各領域の5量で \(|y_{\rm I1b}-(y_{\rm I1a}+2y_{\rm RB1})|\le e_{\rm I1b}+e_{\rm I1a}+2e_{\rm RB1}+8\epsilon(|y_{\rm I1a}|+2|y_{\rm RB1}|)\) |
| I4b の加法性 | 面積 11 列と各領域の \(J_{1,j}\) で \(|y_{\rm I4a}+y_{\rm I4b}-y_{\rm RB1}|\le e_{\rm I4a}+e_{\rm I4b}+e_{\rm RB1}+8\epsilon|y_{\rm RB1}|\)。\(J_0\) は3つの実行でビット単位で同じ |
| Z2 の面積 | \(|A_{\rm Z2}-A_{\rm RB1a}|\le e_{\rm Z2}+e_{\rm RB1a}\)（11 列） |
| RA4b の差 | \(|w_{\rm RA4a}-w_{\rm RA4b}-0.06\ln10|\) が2つの \(w\) の比較幅の和以下 |
| RA1 の乗数 3 | 同じパネルの衛星の各キーで \(|J_3-1.5J_2|\le4\epsilon|1.5J_2|\)（値が 0 の領域では \(J_3=0\) を求める）。面積は層A の結果で、同じ |
| I3、I8、C01、C02 の再生 | ビット単位（配列のバイト列） |
| 成分の関係 | 成功した行で、面積と各領域の `total` と `no_hole`\(+Q_h\cdot\)`hole` の差が、3つの誤差の和と丸め以下。**RA10 には適用しない**（`total` に識別値を注入するため）。収支（\(\sum_DJ_0=S\) など）は `diagnostics.balance` の上界以下。**収支は合格の根拠にしない** |

## 8. 時間の上限

- **CI の増分:** 経過時間で 60 秒以内。Codex の見積もりは、125 点の V で約 43–51 秒と比較の費用（試走の計算の芯 約 21 秒、C04–C08 約 18 秒、N29 約 0.5 秒、V の両版 約 2.2 秒、タグの側の起動など 1–9 秒）。499 点の V では約 53–88 秒で、60 秒を凍結する根拠がない。
- **超えたときの規則:** §5.3 の順（C05、次に C07 を全件だけに回す）。許容幅、期待する状態、全件の対象は変えない。それでも超えれば、旧版を残して再登録する。
- **R-B の全件:** 1コアの経過時間で 2 時間以内。試走の見積もり（CPU 5–10 分）は、この上限に対する見積もりである。
- **追加の点の打ち切り:** 参照の担当が登録の外の点を足すときは、比較を開く前に打ち切る。比較を開いた後に点を足さない。

## 9. 変異の表（必須 16、追加 7）

レビュー担当が、各行の変異を本番のコードに1つずつ入れる。下の表は、登録の内容と入力から決まる判別の条件である。検出の実績は実装の後に測る。対象のキーは表の `mutations` にある。

| ID | 変異 | 登録した行と、落ちるキー・状態 | 区分 |
| --- | --- | --- | --- |
| M01 | 軸の反転 | RA1 の `bin[2]`・`bin[3]` と `outside.high` の `single_loss_total`（表示の衛星が 0、高い側の外が 0.60 になる） | 必須 |
| M02a | 衛星の \(L\) の脱落 | RA1 の `bin[2]`・`bin[3]` の `single_loss_total`（総和 0.30）。\(w\) は変わらない | 必須 |
| M03 | 衛星の乗数 \(L\) を \(1/A\) に置き換える | 同じキー（総和 1.0）。別に計算する \(SLA\) は 0.60 のまま | 必須 |
| M04 | 外側の乗数だけを \(L\) に従わせ、M3a に渡す長さを 20 に固定 | 検出先は RB5: `_call_kernel` の実引数 `weight_length_A`（5 のはず）と `bin[1]`・`bin[2]` の `single_loss_total`。RB1（長さ 20）は変異の影響を受けない正常の対照で、検出先でない | 必須 |
| M05a | `hole_charge` を上流に渡さない（既定の 1 になる） | RB6、RB7 の実引数 `hole_charge` と `bin[1]` の `single_loss_total` | 必須 |
| M06 | 運動量の指定の取り違え | RB1・RB2 の `bin[1]`・`bin[2]`、RB3・RB4 の `bin[0]`・`bin[2]` の `single_loss_total` | 必須 |
| M07 | 一次線のビン幅を密度の換算で余分に掛ける | RA3 の `bin[2].single_loss_hole`（[998,1000]。一次線の第2ビン [999,1002]、幅 3 だけが寄与し、3 倍になる） | 必須 |
| M08a | パネルの半幅の因子の脱落（値の経路） | RA4a の `window_area.term_13a`、`window_area.total` | 必須 |
| M10 | 外を反対の端に足す | RA2 の `outside.low`（0.2）と `outside.high`（0）の `single_loss_total` | 必須 |
| M11a | 外を全体からの残差で求める | RA6 の `outside.low.single_loss_total`（0.6）、RA6b の `outside.high.single_loss_total`（1.2）。§7.2 の置き換えの幅で比べる | 必須 |
| M12 | 内端の δ を両側に数える | RA1 の `bin[5].zero_loss`（0）と `bin[6].zero_loss`（1） | 必須 |
| M13 | 積分後の負の `hole` を切る | RA3 の `bin[2].single_loss_hole`（\(-2/225\) が 0 になる） | 必須 |
| M14a | 初期の 6 枚で打ち切る | RA13 の `window_area.term_13a`、`window_area.total`（\(10^{-30}+10^{-13}|r|\)） | 必須 |
| M14b | 内部の節点の失敗を 0 にする | N29 の終了状態と、stderr の `reason`、`node_kind`、`index`、`loss_eV`、`inner_reason`、終了コード 2、出力先なし。どれかが崩れれば検出 | 必須 |
| M16 | 解像の尺度による事前の分割、上限の検査、節点の分離の検査を迂回する（どれを迂回したかを記録する） | RA8b の終了状態、`exception.required_panels`（1,500,002）、渡した核の評価 0 回。N32b の終了状態、`exception.panel_width_ulp`、`exception.panel_eV`、渡した核の評価 0 回。巨大な配列の確保の失敗は検出に数えない | 必須 |
| M17 | 検証の残差を出力に伝えない | RA12 の `bin[0].single_loss_total` と、細分 1 段以上の条件 | 必須 |
| M05b | `total` を係数 1 で作り直す | RA3 の `window_area.total`（0.165 が 0.18） | 追加 |
| M05c | `total` を正しい \(Q_h\) で作り直す | RA10 の `bin[0].single_loss_total`（3.0 が 0.72） | 追加 |
| M08b | 誤差の経路の半幅の因子の脱落 | N27 の点の誤差の和が \(4\times10^{-4}>3.5\times10^{-4}\) になり、成功が `kernel_error_exceeds_budget` に変わる（予告した状態の変化）。`kernel_point_error` の照合でも落ちる | 追加 |
| M09 | 各セルで \(G\) を中点の値に固定する | RA3 の `bin[2].single_loss_hole`（\(G=(4-T)/3\) と `hole` がともに傾く） | 追加 |
| M11b | 内側のビンを累積の差で作る | RA6 の `bin[1].single_loss_total`（1.2） | 追加 |
| M15 | 片側の感度の左右の取り違え | RA7 の `bin[1].cutoff_sensitivity.plus.total` と `bin[2].cutoff_sensitivity.minus.total`（各 0.40） | 追加 |
| M18 | 点の誤差で重みの絶対値を取らない | RA9 の `bin[0]` の3成分の `kernel_point_error` | 追加 |

**検出の数え方:**

- 同じ条件で正常版が通り、変異を1つだけ有効にして、登録した行の対象のキー、診断、明示した終了状態が落ちたときに、その ID を1回数える。複数の行やキーが落ちても足さない。
- 判定は全変異で同じである: **表の `detected_by` のキーのいずれかを破れば検出**（表の `state_change`）。状態の変化は、「成功を返した場合」に限らず、対象とした拒否、その属性、評価の前の停止のどれを破ったかで見る。予告した計算上の状態の変化（N27 など）と、試験器の例外や無関係なクラッシュを区別し、後者は数えない。
- 表の `unaffected_controls` は、変異の後も正常版と一致するはずの正常の対照である（M04 の RB1）。
- skip、比較 0 件、変異が入っていない場合、参照なし、正常版の失敗、登録の外の検査だけの失敗は、登録先での検出に数えない。
- 報告は、必須 \(k/16\)、追加 \(k/7\) とし、生存、未実行、判定不能を明示する。
- 計画の \(L\) の二重の変異（M02b）は、M02a と同じ検出先の反対向きの誤りなので外した（U17 で承認）。

## 10. 凍結の前の確認

**下書きの担当が確かめたこと（2026-11-26）:**

- 材料のファイルの SHA-256（Al、Si）は、`benchmarks/inputs/` の現在のバイト列と一致する。
- §1.3 の16進は、十進の式を binary64 で計算した値と一致する（N08、N09 は設計 §6.5 の16進と同じ）。
- §1.4 の全件の GA と GC は、標準ライブラリだけで作った配列の SHA-256 が M3a の事前登録の値と一致した。CI の2つの配列は同じ作り方で作った。
- 表の生成は、100 行、正規化した入口の一意性（説明だけの欄を除き、API の入口は実効値で比べる。重複として拒否されるべき反例の自己検査を含む）、族ごとの数、CI の 72 行、理由の語彙、再利用の向きと循環のないこと、予定のキーの数（R-A 922、CE 66、R-B 1048）、変異の ID の集合と対象のキーが台帳にあることを、例外で検査する。`--check` はバイト単位で一致し、ハッシュの種を変えても同じだった。

**Codex の1回目の確認（2026-11-26、第1稿に対して）:** 判定は「修正後に凍結」。blocker 2、major 5、minor 2。全文は `work/consult/m4a/prereg/codex_prereg_review.md`。Codex は、V01–V08 の入力の状態と時間、N08・N23・N28 の上流の状態、表の再現（第1稿の SHA-256 `327ea2a1…`）を確かめた。指摘はすべて Codex の案のとおりに採った。設計 第6版と U1–U19 に反するために採らなかった案はない。

| 指摘 | 対応 | 場所 |
| --- | --- | --- |
| 1 [blocker] 故障を注入した層B に核の合格条件を課している | RA9・RA10 の成功を「固定パネルで層B が返る」と定め、駆動の成功の判定を外した。RA10 は成分の恒等式の適用外。RA9 は添字 2（左端が 0）の節点の \(|\omega_2|Le\) を `kernel_point_error` と比べる | §5.1、§5.2、§7.3、§7.4、表の `expect` |
| 2 [blocker] 比較のキーの分母と判定の方法が決まらない | 各行に `check_plan`（キー、検査の種類、オラクル）を付け、台帳を展開した。必須の数値の比較は 2036 キー（R-A 922、CE 66、R-B 1048）。I1a・I4a は照合済みに数えない。逆数は2通りの検査で、§7.3 を重ねない。\(r=0\) の順と、`no_reference` などを台帳に残す規則を書いた | §5.1、§6.3、§7.1–§7.3、表 |
| 3 [major] 100 の ID と独立な入口の数が違う | V07 の主の実行を P1 の深さ平均だけにし、残りの2つを補助の対照にした。V08 は `experimental` を求める。`reuse` を成果物まで書き、C04–C08 の再利用を外した。`execution_order` を表に入れた | §5.1、§5.2、表 |
| 4 [major] 返値の読み方と失敗時の観測先が足りない | 返値の欄の名前と形、`error_budget.outputs[key]` と同率の順、例外の `diagnostics` のスナップショットと属性の表、CSV の見出しと順、`run.json` のパス、`spectrum.json` の形、精度の優先順、CONFIG の `integration` の明示、`resolution_basis` の `null`、継承した警告の規則、`capabilities` に足す識別子を書いた | §2、§4.1–§4.5 |
| 5 [major] 解像の規則を直接に観測できる条件 | 分離の3つの場合と、208 ulp を閾値にしないこと。根の変数、確認、相殺、未確認の根の保持、`end: both`。CE の族と端の指定、CE3 の条件つきの警告、計算の `area_epsabs` と比較の床の併用。領域ごとの初期のセル。点の数は1つのメッシュの数で、累積の上限でないこと、`kernel_calls` の欄 | §2、§3、§4.4、§5.2 |
| 6 [major] CI 60 秒の判断 | V を CI で 125 点、全件で 499 点にし、配列の SHA-256 を別々に凍結した。60 秒を超えたら C05、次に C07 を全件だけに回す。上限は 60 秒のまま | §1.4、§5.3、§8、表の `v_profiles`、`ci_reduction_order` |
| 7 [major] 変異の検出の経路 | 変異ごとに、登録した行、落ちるキー、状態の変化を書き、表の `mutations` に対象のキーを入れた。検出の数え方を改めた | §9、表 |
| 8 [minor] 表の再現性と集計の意味 | 生成の検査に、正規化した入口の一意性、キーの数、変異の集合と対象のキー、再利用の循環を足した。RA8a の必要数 41 を期待に入れた（RA13 の 6 も） | §5.2、§10、生成のスクリプト |
| 9 [minor] 未決の表の軽微な項目 | \(u_r\) と丸めの単位の記号を分けた。節点の添字は左端を 0、解析の層は \(E_0\) を持たない、N11 は2種類の一次線の組、CE の候補は族・端・別々の根で照合、L は実行して記録した件数、とした | 冒頭、§3.1、§4.7、§5 |
| A1–A6、B1、B2 | Codex の決定案のとおり。A1・A2 は §7.2、A3 は §4.6、A4 は §4.7、A5・A6 は §5.2・§5.3・§1.4、B1 は §8、B2 は §7.2 に書いた | — |

**Codex の2回目の確認（2026-11-26、第2稿に対して）:** 判定は「修正後に凍結」。凍結前の修正は major 5、minor 2。全文は `work/consult/m4a/prereg/codex_prereg_review2.md`。1回目の指摘は、1・3・6・9 と A6・B1・B2 が解消、残りが一部と判定された。Codex は、表の再現（第2稿の SHA-256 `fb9f85df…`、ハッシュの種 1 と 77）、キーの数（全 2,982、必須の数値 2,036）、4つの配列の SHA-256、タグのコミットを確かめた。§4 の具体の名前（R1）と `error_part`・RA1 の幅（R4）は採用された。修正はすべて Codex の案のとおりに採った。

| 修正 | 対応 | 場所 |
| --- | --- | --- |
| 1 [major] C01 の再生と I3・Z2 の観測先 | C01 の観測先を `spectrum_definition.config_loss_grid_used` に統一した。設定を面積 6 個・スペクトル 8 個に分けた。N29 の CONFIG と argv を表と本文に書いた。`diagnostics.primary_support`（4つの値）を置き、マニフェストの `primary.support` に写す。C01・C02 は `results` と CSV の全数値の列を API と照合する | §4.1、§4.4、§5.2、表 |
| 2 [major] 層B の値と誤差の取得 | `fold` の戻り値を `folded.outputs[key]`（`value`、誤差の5項、`estimate`）にした。N26・N27 の `initial_edges=[2,4,5]` を明記した。N28・N29 は `outputs` が空、`max_ratio` が `None`。解析の行の「評価 0 回」は、渡した `kernel(T, side)` の呼び出しで数える | §4.4、§4.7、§5.2、表 |
| 3 [major] ゼロ参照と完了表の欠値 | \(r=0\) の比較幅を \(d_{\rm abs}\) とした。逆数の (i) は同じ行の `window_area.total` の参照を使う。§4.6 は、全条件が一致した記録で全必須キーが `pass` のときだけ `matched`、ほかは `not_matched` と確定した（`no_reference` の明記を含む）。`_load_comparison_table()` の戻り値の形を決め、C08 に `no_reference` を明記した表の行内の検査を足した | §4.6、§6.3、§7.2、表 |
| 4 [major] V のタグ比較に混ざる正当な差 | V05・V06 の stdout の `output_dir` は各版で確かめてから外す。V08 は `version` を各版の `__version__` と照合し、接頭部・既存のブロック・追加の2名を別に検査する。V07 の Si の補助の対照を `np.r_[GA_full[:10], 1.12]` と書いた | §5.2、表 |
| 5 [major] 変異の指定 | M04 の検出先を RB5 だけにし、RB1 を正常の対照（`unaffected_controls`）にした。M14b に stderr の添字・損失など、M16 に `exception.required_panels`・`panel_width_ulp`・`panel_eV` を足した。状態の変化の判定を書いた | §9、表 |
| 6 [minor] RA6・RA6b の第1ビンの置き換えの強度 | \(2^{60}\) に直し、駆動の全強度 \(2^{60}+3\) と書き分けた | §7.2、表 |
| 7 [minor] 入口の一意性のキー | 説明だけの欄（`input_template`、`read_from`）を除いて検査する | §5.1、生成のスクリプト |
| R1 の追加の提案 | `dependent_model_ids` を3つのキーの辞書に、継承した警告に `occurrences` を、`batch_sizes` の定義を採った | §4.4、§4.5 |
| R4 | RA1 の式を \(|J_3-1.5J_2|\le4\epsilon|1.5J_2|\) にした | §5.2、§7.4 |

**Codex の3回目の確認（2026-11-26、第3稿に対して）:** 判定は「修正後に凍結」。2回目の修正 1・2・4・6 は解消、3・5・7 は一部と判定された。残りは「凍結前に揃える箇所」の3点（major 1、minor 2）。全文は `work/consult/m4a/prereg/codex_prereg_review3.md`。Codex は、表の再現（第3稿の SHA-256 `87532028…`）、キーの数（必須の数値 2,036）、CI 72 行と全件だけ 28 行を確かめ、第3稿が置いた名前（`comparison_table_shape` の最上位の3つ、`test_table_C08_no_reference`、`kernel_evaluation_counters`、`replacement_intensity_RA6`、`defaults.area_settings`・`spectrum_settings`）を採用した。3点はすべて Codex の案のとおりに採った。

| 揃える箇所 | 対応 | 場所 |
| --- | --- | --- |
| 1 [major] M04 と M14b の判定文 | §6.2 の幅の超過は RB5 だけに課し、RB1 には「変異の後も正常版と一致する」条件を書いた（表の `unaffected_controls` の `condition`）。全変異の `state_change` を「`detected_by` のいずれかを破れば検出」に揃えた（M14b の `stderr.inner_reason` の抜けがなくなる）。個別の注記は `note` に分けた | §6.2、§9、表の `mutations` |
| 2 [minor] C02 の比較の条件 | 正規化した `conditions` の JSON の例を1件固定した（キーの名前、入れ子、配列の順）。面積の入口とビンの一次線の違いも書いた | §4.6、表の `conditions_example_C02`、`conditions_rules` |
| 3 [minor] 入口の重複の検査 | API の入口は、既定の物理の引数と数値の設定を展開した実効値で比べる。CLI の差し替える表と C06 の生の入力の差はキーに残す。自己検査に、重複として拒否されるべき反例を入れた: RB5 の長さを明示の 20 Å にした表（RB1 と重複）、N36 を `max_cells=88` と明示の `epsrel=1e-6` にした表（N30 と重複）、C08 の表を C02 と同じにした表（C02 と重複）。C06 の表だけを揃えても重複にならないことも確かめる | §5.1、生成のスクリプト |

**4回目の確認はしない。**統括担当が、この3点の反映と `--check` を確かめて凍結する。表のスキーマは v4 になった。100 行の入力、キーの数、公開の成功の条件、時間の上限は第3稿から変えていない。

**統括担当の凍結の記録（2026-11-26）:** 第4稿の3点の反映（M04 の検出先を RB5 にし RB1 を正常の対照にしたこと、全変異の `state_change` の一本化、C02 の条件の例、実効値での重複の検査と自己検査の反例）と、`--check` がハッシュの種を変えてもバイト単位で一致することを確かめて凍結した。表の JSON の SHA-256 は `250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408`（スキーマ v4、100 行、CI 72 行、必須の数値の比較 2,036 キー）。凍結は、この文書と表を main に入れたコミット（題「M4a: design decided and checks pre-registered」）である。未決の一覧の残り（CI の実時間、タグとの一致、実装の後にしか観測できない状態と変異の検出の実績）は、実装の後に記録する。

## 付録　参照値の凍結

（参照の担当が、実装の結果と比べる前に記入する。）

### R-A

- `benchmarks/numerics/m4a/references_ra_m4a.json` の SHA-256: `＿＿＿＿`
- 生成のスクリプトとその SHA-256、コミット、環境、所要時間:
- キーの数（行ごと。台帳の数との照合）、\(u_r\)/幅 の最大:
- RA13 の凍結前の確認（正常版と初期打切り版の判定）:
- 既知の誤答の記録:
- 「値なし」のキーと、当てはめた条件（§6.3 の番号）:

### R-B と CE

- `benchmarks/numerics/m4a/references_rb_m4a.json` の SHA-256: `＿＿＿＿`
- 生成のスクリプトとその SHA-256、コミット、環境、所要時間（1コアの経過時間と CPU）:
- キーの数（行ごと。R-B 1048、CE 66 との照合）、\(f_r\) の最大、\(u_r\) の内訳の最大:
- 識別のビンの確認（2指定の差と比較幅の和の比）、M04・M05a の判別の確認:
- 候補の根を別に計算したことの記録:
- 「値なし」のキーと当てはめた条件、参照の経路の失敗:
