# xpsloss

公開理論に基づくXPSの誘電損失計算を、Python APIとCLIで再現するためのパッケージです。0.3.0 の対象は M0、M1、M2-1、硬いギャップのある材料の固定深さ（M2-2）、外から与える一定の長さでの法線出射の深さ平均（M3a）と、法線出射の単一損失スペクトルと有限の損失窓の面積（M4a）、その重みの長さと乗数を分ける規約（M4a-2）です（v0.1.0 は M0、M1、M2-1。0.2.0 で M2-2 と M3a を加えた）。開発版（0.4.0.dev0、未公開）には、法線出射の多重損失の近似の模型（M4b-1、下記）が入っています。[English README](README.md)

**現時点で動くのは、次の6つです（6つ目は 0.3.0 から）。開発版（0.4.0.dev0）には7つ目と8つ目があります。**
- **振動子で表した損失関数**
- **その因果的な複素誘電応答**
- **均質・等方な無限媒質のバルク非弾性散乱**
- **法線出射・固定深さでの表面・正孔・真空中の項（Simonsen 1997 の公表式。ギャップのない材料と、0.2.0 からは硬いギャップのある材料）**
- **法線出射の深さ平均（原著の式(18)、外から与える一定の長さ）**
- **法線出射の単一損失スペクトルと有限の損失窓の面積（原著の式(19)。重みの長さと乗数を選べる）**
- **（開発版 0.4.0.dev0）法線出射の多重損失の近似（最初の損失は全核、2回目以降は 13a の項だけ。次数の上限 `max_order` は必須）**
- **（開発版 0.4.0.dev0）その多重損失の束縛エネルギーの入出力と測定の記録（入口と出口の変換だけ。光子エネルギーは必須。パッケージは軸をずらさない）**

**斜出射、軌道ごとの長さや自己無撞着な長さでの深さ平均は未実装です。多重散乱は、開発版に上の近似の模型が1つあるだけで、物理の妥当性は確かめていません。QUEELS-XPSの完全再現版ではありません。分散のあるAlでは、QUEELS-XPS v2.1 との比較で、両方の運動量指定が事前の基準を満たしませんでした。**

このパッケージはQUEELS/QUASESの公式製品ではなく、その実行ファイル・DLL・ソースコードを利用していません。公開されている数式から新規に実装しています。

## 現在の範囲

| 機能 | 状態 |
| --- | --- |
| JSON入力の検証 | 実装済み |
| Drude-Lindhard型振動子からELFを計算 | 実装済み |
| 運動量分散を含むバルクDIIMFP | 実装済み |
| CSV、JSONによる出力、計算条件の記録 | 実装済み |
| 解析解・和則・収束との照合 | テストあり |
| 振動子ELFからの因果的な複素誘電応答（ギャップ時はKKで実部を再構成） | 実装済み・テストあり・別担当がレビュー |
| 表形式ELFからのKK再構成（参照・診断用） | 実装済み・テストあり・別担当がレビュー |
| 光学和則・一次モーメントの診断（入力は再規格化しない） | 実装済み・テストあり |
| 文献由来の固定入力（ELFデータ集のSi・Au・Al）と来歴の記録 | 同梱。値の物理的妥当性は未検証 |
| 法線出射・固定深さの表面・真空中の損失・静止正孔・干渉（ギャップなし、[Simonsen 1997 式(13)–(17)](docs/M1_EQUATION_MAP.md)） | 実装済み・テストあり・別担当のエージェントがレビュー。原著の図3・5・8は、別担当のエージェントが原著の式(20)を使った独立なコードで再現を試みた。比較点84点のうち、読取り誤差の範囲に入ったのは36点、超えたのは18点、幅の仮定に依存したのは30点（原著の Table I に基づく一律の幅の規則による。[比較点ごとの結果](benchmarks/paper_figures/RESULTS.md)） |
| 同上の QUEELS-XPS v2.1 との照合 | 分散0の合成入力では、成分ごとの最大絶対差は参照 `total` の最大値の 0.1–0.4 %。分散ありAlの結果と合成入力の追試は下記。 |
| 法線出射・固定深さ、硬いギャップのある材料（0.2.0、M2-2） | 実装済み・テストあり・別担当がレビュー。ELF データ集の Si と合成入力で、実装担当とは別の担当がパッケージを使わずに多倍長で計算した参照値194点と、事前登録の許容幅で一致した。ギャップ下は、減衰のない表面モードがないことを証明できる損失だけを計算する |
| 斜出射 | 未実装 |
| ギャップ下の減衰のない表面モード（ε = −1）の寄与、ギャップ端の広がり | 未実装（該当する損失は拒否） |
| 法線出射の深さ平均、外から与える一定の長さ（0.2.0、M3a） | 実装済み・テストあり・別担当がレビュー。実装担当とは別の担当がパッケージを使わずに多倍長で計算した参照値498キーと、事前登録の許容幅で一致した。事前登録の全2041実行が、登録した状態（成功、拒否、未収束）になった。公開の固定深さを数値で深さ求積した経路とは、値のある630点で一致し、26点は経路の値がなく未完了。`k_parallel` の独立経路 I3 との照合（V9）は計算量のため行っていない |
| 法線出射の単一損失スペクトルと有限の損失窓の面積、外から与える一定の長さ（0.3.0、M4a。重みの長さと乗数の規約は M4a-2） | 実装済み・テストあり・別担当がレビュー。原著の式(19)のゼロ損失線と単一損失の衛星を、損失窓の中だけで計算する。事前登録の 100 行が登録した状態になり、実装担当とは別の2つの担当が作った参照の 2041 キー（解析の参照 927、ELF データ集の Al・Si の通し経路 1048、動く端の合成入力 66）と、事前登録の許容幅で一致した。成功は推定誤差が予算の内という意味で、証明された上界ではない。実測スペクトルとは比べていない |
| 法線出射の多重損失の近似 `first_loss_13a_chain_v0`（開発版 0.4.0.dev0、M4b-1。計算の軸は運動エネルギー） | 実装済み・テストあり・別担当がレビュー。最初の損失は M4a の全核、2回目以降は 13a の項だけを繰り返す。事前登録の 109 行が登録した状態になり、実装担当とは別の担当が作った参照（解析の参照 3,572 キー、ELF データ集の Al・Si の参照 6,096 キー）と、事前登録の許容幅で一致した。確かめたのは、与えた核と模型の定義に対する畳み込みの数値で、模型の物理の妥当性は確かめていない。実測スペクトルとは比べていない |
| 軌道ごとの長さ・自己無撞着な長さでの深さ平均 | 未実装 |
| 上の多重損失の束縛エネルギーの入出力と測定の記録（開発版 0.4.0.dev0、M4b-2） | 実装済み・テストあり・別担当がレビュー。フェルミ準位基準の束縛エネルギーを、入口と出口で E = C − B と換算するだけで、物理は足さない。光子エネルギーは必須で、パッケージは軸をずらさない。事前登録の 82 行が登録した状態になり、実装担当とは別の担当が有理数で作った参照（1,771 キー）を含む台帳の 25,948 キーがすべて合格した。模型の物理の妥当性は確かめていない。実測スペクトルとは比べていない |
| 実測スペクトルとの比較（M4b-3） | 未実装 |
| 上の近似のほかの多重散乱の模型、窓の外の損失、内因性の損失の分離、装置関数、背景の差引き | 未実装 |
| Quanty接続、CoO物性値、フィッティング | 未実装 |
| 元のQUEELS-XPSとの数値一致 | 主張しない（上記） |

分散ありAlでは、電子エネルギー1000 eV、法線出射、深さ9.967 Å、損失1.0–99.8 eVで両指定が不合格でした。参照 `total` の最大値に対する `total` の最大絶対差は `k_total` で6.6 %、`k_parallel` で7.5 %でした。合成3入力の追試では、各実行で深さが10 Åに最も近い列の損失5.0–99.8 eVを比べました。`k_total` と候補式を使うと、全3成分の差がその列・損失区間の参照 `total` の最大値の0.3 %以内に入りました。本体の公表式では各入力で少なくとも1成分が基準を超えました。残る約10⁻³の差の原因と参照ソフトの内部処理は未確認です。[Alの比較](docs/M2_FIRST_CASE.md)・[合成入力の追試](docs/QUEELS_DISPERSION_TEST.md)

ELFは Im[-1/epsilon]、DIIMFPは単位経路長・単位損失エネルギー当たりの非弾性散乱確率密度です。ELF単体を、測定されるXPSの損失確率と同一視しないでください。

## Windows PowerShellでの開始

Python 3.10以上を使います。展開した `xpsloss` のルートフォルダーを開きます。

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m xpsloss capabilities --json
.\.venv\Scripts\python.exe -m xpsloss validate examples/synthetic_bulk.json
.\.venv\Scripts\python.exe -m xpsloss elf examples/synthetic_bulk.json --out runs/elf-demo
.\.venv\Scripts\python.exe -m xpsloss dielectric examples/synthetic_bulk.json --q-inv-A 0.5 --out runs/dielectric-demo
.\.venv\Scripts\python.exe -m xpsloss bulk examples/synthetic_bulk.json --out runs/bulk-demo
.\.venv\Scripts\python.exe -m xpsloss xps-fixed-depth examples/synthetic_bulk.json --depth-A 10 --epsilon-momentum k_total --out runs/xps-fixed-depth-demo
.\.venv\Scripts\python.exe -m xpsloss xps-depth-averaged examples/synthetic_bulk.json --weight-length-A 20 --epsilon-momentum k_total --out runs/xps-depth-averaged-demo
.\.venv\Scripts\python.exe tests\run_suite.py
```

`xps-fixed-depth` は、表面法線からの出射角0°、生成深さ `--depth-A` で、正孔ありの有効断面積と、その「正孔なし」「正孔（単位電荷）」成分、項ごとの値を出力します。`--epsilon-momentum` は必須です。`k_total`（全運動量で ε を評価）を選ぶと、深い極限でバルクの式に一致します。`k_parallel` は原著の閉形式の導出どおりの評価です。

2つの指定は、確かめ方が違います。
- `k_parallel`: 独立な数値経路 I3 で、分散のない入力と分散のある入力の全8項を照合しました。
- `k_total`: 13a のバルクの恒等式、深い極限、同じ閉形式の別の実装との照合で確かめました。分散のある入力で、有限深さの全項を元の電荷・電位の式から独立に検証したものはありません。
- 分散のない入力では、2つの指定は一致します。

計算は原著の仮定に従います。t = 0 に電子と正孔が同じ点に突然生じ、正孔は静止して寿命は無限です。電子は直線の軌道を、初期運動エネルギーで決めた一定の速度（非相対論）で進みます。損失は E0 に比べて小さいとし、表面は鏡面反射の模型で扱います。

`run.json` の `applicability` には、最大損失と E0 の比、モデルの仮定、照合範囲とその根拠が入ります。E0 の範囲は 1000–4000 eV、比の上限は 0.1 です。E0 の両端で、T/E0 が 0.1 までの離散点を使い、`k_parallel` の全8項を独立経路と照合しました。深さ、材料、線幅、分散、正孔電荷、損失格子の下端は範囲判定の対象外です。どちらかの範囲を外れた計算が成功すると、`warnings` に外れた量を記録します。範囲内でも、任意の入力での精度や物理的な妥当性は未検証です。

硬いギャップのある材料（M2-2、0.2.0 から）では、段差の ELF とクラマース・クローニッヒ変換による実部を使い、広がりは加えません。ギャップ端のごく狭い除外幅の損失は拒否します。ギャップより下では媒質に損失がなく、運動量の範囲に減衰のない表面モード（ε = −1）がないことを外向きに丸めた解析的な上界で証明できた損失だけを計算し、それ以外は `FixedDepthGapDomainError` で拒否します。そこでは正孔電荷1の `total` が 0 になり、`no_hole` と `hole` は吸収ではない反応的な成分です（警告に記します）。振動子ごとの単調性の条件も受理の条件です。ギャップありで照合した範囲は、ELF データ集の Si、1000 eV、深さ 2–200 Å、損失 0.2–99.8 eV で、それ以外の入力には警告を出します。[M2-2 の設計](docs/M2_2_DESIGN.md)・[事前登録](docs/M2_2_PREREG.md)

次の2つは 0.3.0 からあります。`examples/single_loss_spectrum.json` は合成の一次線（1000 eV の δ 線1本）と出力の端です。

```powershell
.\.venv\Scripts\python.exe -m xpsloss xps-loss-window-area examples/synthetic_bulk.json --loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total --out runs/loss-window-area-demo
.\.venv\Scripts\python.exe -m xpsloss xps-single-loss examples/synthetic_bulk.json --spectrum examples/single_loss_spectrum.json --loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total --out runs/single-loss-demo
```

深さ平均（M3a、0.2.0 から）は、`effective_cross_section_depth_averaged` と `xps-depth-averaged` で、法線出射の固定深さの11出力を、原著の式(18)の単一散乱の重み x exp(−x/L)/L² で深さ方向に平均します。重みの長さ L（`weight_length_A`、`--weight-length-A`）は必須で、呼び出す側が与え、全経路で一定です。出所は `--weight-length-source` に書けます。深さの積分は運動量の被積分関数の中で解析的に閉じるので、深さの格子はありません。11出力は固体内の経路長当たりの Å⁻¹ eV⁻¹ で、真空中の損失を含み、`total = no_hole + hole_charge * hole` が成り立ちます。原著の留保（式(18)は REELS の重みの流用で XPS での妥当性は検討を要すること、結果は最初の1回の非弾性散乱だけを表し intrinsic の部分を含むこと）は `warnings` に書きます。軌道ごとの長さや自己無撞着な長さ、弾性散乱、一様でない深さ分布、損失窓の面積、スペクトルは扱いません。未収束は `DepthAverageConvergenceError`（`reason` つき）です。事前登録の検査で使ったのは、ELF データ集の Al・Si と合成入力、ほぼすべて 1000 eV、長さ 5、20、80 Å と、極限の確認の 10⁻⁶–10⁶ Å、両方の運動量の指定です。`run.json` は照合した条件を記しますが、実行ごとの範囲の判定はしません。[M3a の設計](docs/M3A_DESIGN.md)・[事前登録](docs/M3A_PREREG.md)

単一損失のスペクトルと損失窓の面積（M4a、0.3.0 から）は、`loss_window_area`、`single_loss_spectrum` と、`xps-loss-window-area`、`xps-single-loss` で計算します。原著の式(19)を、法線出射、比例定数 1、外から与える一定の長さ L で使い、ゼロ損失線 J0 = F と単一損失の衛星 J1 = L∫F(E+T)K_sc(E0,T;L)dT を返します。核は深さ平均（M3a）の出力で、すべての一次線のエネルギーに対して設定の `electron_energy_eV`（E0）に固定します。既定では、L を式(18)の重みと式(19)の乗数の両方に使います。両者は別々に選べます（M4a-2、下記）。

- **損失窓は必須です。** `--loss-window-eV a b`（`loss_window_eV`、0 < a < b < E0）に既定値はありません。原著は式(19)の積分が発散するとして下限 0.05 eV から積分しており（p. 1618）、結果は下限に依存しえます。窓の外の損失は計算せず、収支にも入れません。下限への感度（`cutoff_sensitivity`）を一緒に返します。
- **面積:** `window_area.csv` に、11列（8項、`no_hole`、`hole`、`total`）の窓の面積（Å⁻¹）、推定誤差、予算、下限への感度（無次元）が入ります。`run.json` には `single_loss_weight`（L×`total` の面積、無次元。確率ではなく 1 を超えうる）と、面積の逆数とその区間（Å。平均自由行程とは呼びません）が入ります。
- **一次線の形式:** `--spectrum` の JSON に、`schema`（`xpsloss-single-loss-spectrum-v0`）、`axis`（`kinetic_energy_eV`）、`primary`、`output_edges_eV` を書きます。`primary` は、δ 線（`kind: "delta"`、`energies_eV`、`intensities`）か、ビン（`kind: "binned"`、`edges_eV`、ビンごとの積分強度 `intensities`）のどちらか一方です。強度は非負で、位置は非負の運動エネルギーです。出力の端は配列か `{"min", "max", "step"}` で与えます。一次線は内因性の損失の衛星を含まないものとして扱います。
- **ビン積分と単位:** `single_loss_spectrum.csv` は、出力のビンごと（左閉右開、最後のビンだけ右端を含む）と両側の外（`low`、`high`）について、ゼロ損失、衛星の3成分（`total`、`no_hole`、`hole`）、その和 `spectrum` のビン積分と推定誤差を、一次線の強度の単位で返します。負の値は切りません。密度（強度の単位/eV）は `run.json` の `*_per_eV` にあります。設定ファイルの `loss_grid_eV` は検証しますが使いません。
- **成功の意味:** 成功は、面積と各ビンの推定誤差が予算の内という意味で、証明された誤差の上界ではありません。解像の根拠の強さは `diagnostics.resolution_basis` に記録します。未収束は `SingleLossConvergenceError`、入力の誤りは `SingleLossInputError`（どちらも `reason` つき）で、CLI は終了コード2で止まり、出力先を作りません。
- **重みの長さと乗数（M4a-2）:** 式(18)の重みの長さ L_w と式(19)の乗数 m を別々に選べます。J1 = m∫F(E+T)K_sc(E0,T;L_w)dT で、単一損失の重みは m と `total` の窓の面積の積です。重みの長さは、外から与える値（`--weight-length-A`）か、`--weight-length-inverse-13a-window-area` による L_w = 1/(A_13a + c) です。A_13a は 13a の項の窓の面積（まず種の長さ `--weight-length-seed-A`、既定 20 Å で計算します）、c ≥ 0 は呼び出す側が与えてよい追加の逆数の長さ（`--extra-inverse-length-per-A`）です。この長さは窓に依存し、平均自由行程ではありません。乗数は L_w（既定）、外から与える値（`--multiplier-A`）、`total` の窓の面積の逆数（`--multiplier-inverse-window-area`。単一損失の重みは定義により 1）のどれかです。API のキーワードは `weight_length_kind`、`extra_inverse_length_per_A`、`weight_length_seed_A`、`multiplier`、`multiplier_A` です。M4a の模型の識別子、数値、照合の表示を保つのは、既定の選び方と、L_w に等しい外部の乗数だけです。ほかの組み合わせは `simonsen1997_single_loss_spectrum_conventions_v0` を返し、照合は常に `not_matched` で、警告が増えます。重みの長さと乗数が違う組み合わせは、式(19)の導出の外にあります。窓の下限への感度は L_w と m を固定した値です。乗数に区間があるとき、返す誤差はその分を含むので、予算のおよそ2倍までになりえます。重みの長さの不確かさが核を通して効く分は含みません。[M4a-2 の設計](docs/M4A2_DESIGN.md)・[事前登録](docs/M4A2_PREREG.md)
- **照合の表示:** `run.json` の `comparison` の `matched` は、入力の条件が、配布する完了した比較の表の登録した条件（ELF データ集の Al・Si、登録した一次線・出力の端・窓・長さ・精度）と数値として完全に一致し、その必須キーがすべて合格だったことだけを意味します。核の条件だけが合えば `kernel_conditions_checked`、ほかは `not_matched` です。上の合成の例は登録した条件の外なので、`matched` にならず、警告が出ます。`physical_validity` は常に `not_established` です。
- **範囲の外:** 多重散乱、斜出射、軌道ごとの長さや自己無撞着な長さ、弾性散乱、窓の外の損失、内因性の損失の分離、装置関数、背景の差引き、フィットは扱いません。活動線のある材料で、ギャップ端をまたぐ窓とギャップ下の窓は拒否します。QUEELS-XPS の完全な実装ではありません。

[数式と単位](docs/EQUATIONS.md)・[M4a の設計](docs/M4A_DESIGN.md)・[事前登録](docs/M4A_PREREG.md)・[検証記録](docs/VALIDATION.md)

### 多重損失（開発版 0.4.0.dev0、M4b-1）

**公開した 0.3.0 には入っていません。** 開発版のソースツリーで使えます。`multiple_loss_spectrum`、`multiple_loss_order_areas` と `xps-multiple-loss` は、単一損失（M4a）に、同じ損失窓の中の損失を繰り返した次数 n ≥ 2 の寄与を足します。模型は `first_loss_13a_chain_v0` の1つだけで、近似です。法線出射だけを扱います。計算の軸は運動エネルギーです。束縛エネルギーの軸での入出力は、下の M4b-2 の節にあります。

```powershell
.\.venv\Scripts\python.exe -m xpsloss xps-multiple-loss examples/synthetic_bulk.json --spectrum examples/single_loss_spectrum.json --loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total --model first_loss_13a_chain_v0 --max-order 3 --epsrel 1e-7 --out runs/multiple-loss-demo
```

```python
import numpy as np
from xpsloss import DeltaPrimary, Material, Oscillator, multiple_loss_order_areas, multiple_loss_spectrum

material = Material("synthetic_single_oscillator_not_real_material", 0.0, (Oscillator(15.0, 180.0, 2.0, 0.5),))
primary = DeltaPrimary([1000.0], [1.0])
edges = np.arange(940.0, 1001.5, 1.0)
common = dict(loss_window_eV=(1.0, 50.0), weight_length_A=20.0, epsilon_momentum="k_total",
              model="first_loss_13a_chain_v0", max_order=3)

result = multiple_loss_spectrum(material, 1000.0, primary, edges, epsrel=1e-7, **common)
print(result.truncation["state"], result.truncation["support_order"])
for order in result.orders:
    print(order.order, order.area["total"], order.area_error["total"])

areas = multiple_loss_order_areas(material, 1000.0, epsrel=1e-7, **common)
print([order.area["total"] for order in areas.orders])
```

- **模型の近似:** 最初の損失には M4a・M4a-2 の全核（表面・真空・正孔・干渉の項を含む 11 列）を使い、2回目以降の損失には 13a の項だけを使います。正孔・干渉・表面・真空の項は、軌道に高々1回しか入りません。根拠は原著の留保です。式(18)の核は最初の非弾性散乱だけを表し、内因性の部分を含むので、繰り返しの畳み込みで全スペクトルを作れません（Simonsen ら 1997, p. 1615）。全核をそのまま繰り返す模型は公開していません。軌道ごとの深さの相関は落としています。
- **乗数と核:** 繰り返しのたびに、1次と同じ乗数 m を掛けます（n 次は mⁿ に比例）。重みの長さと乗数の選び方は M4a-2 と同じです。核は、全次数で設定の `electron_energy_eV`（E0）に固定します。次数ごとにエネルギーを下げません。
- **`max_order` は必須です。** `--max-order`（`max_order`）と `--model`（`model`）に既定値はありません。`max_order` は 1 以上の整数か `display_support`（表示のビンに届く次数まで全部）です。面積だけの `multiple_loss_order_areas` は整数だけを受け付けます。
- **打ち切りの状態と上界:** `truncation.state` は、表示のビンに届く次数をすべて計算したら `complete_on_display`、そうでなければ `truncated` です（`truncated` も成功として返し、警告を出します）。落とした次数の上界 `truncation.neglected_at_or_above_region_bound` は、ビンごとの値ではなく、そのビンとそれより高い運動エネルギーの側の全部に落ちる量の累積の上界で、推定の誤差を入れた上界です（数学的な保証ではありません）。`--neglect-tolerance-fraction`（既定値なし。一次線の絶対の強度に対する比）を与えると、表示のビンすべての上界が許容以下なら `bounded_on_display`、超えれば失敗になります。上の例は、表示に届く次数が 60 で `max_order` が 3 なので `truncated` です。
- **出力:** `multiple_loss_spectrum.csv` は、ビンと両側の外ごとに、ゼロ損失、次数ごとの3成分（`order<n>_total`、`order<n>_no_hole`、`order<n>_hole`）、その和 `spectrum` と推定誤差を返します。`order_areas.csv` は、次数ごとの 11 列の面積（無次元。確率ではありません）と、窓の下限への感度です。1次の欄は M4a・M4a-2 の結果とビット単位で同じです。高次の誤差は、1次から継承した分、途中の次数から伝わった分、その次数の網の分、乗数の区間の分に分けて返します。負の値は切りません。
- **精度の設定:** 既定の `epsrel`（1e-6、M4a と同じ）では、`max_order` が 2 以上の典型の実行（ELF データ集の Al、δ の一次線、1 eV のビンなど）が、2次で `inherited_error_exceeds_order_budget` で止まります。1次の推定誤差を高次に引き継いだ分が、高次の予算を超えるためです。**`epsrel=1e-7`（CLI は `--epsrel 1e-7`）を指定すると通ります。** 上の例もこの値を指定しています。この設定での Al の例は 34 秒で、同じ引数の M4a は 16 秒でした。
- **成功の意味と失敗:** 成功は、推定誤差が予算の内という意味で、証明された誤差の上界ではありません。未収束は `MultipleLossConvergenceError`、入力の誤りは `MultipleLossInputError`（どちらも `reason` つき）で、CLI は終了コード2で止まり、出力先を作りません。M4a・M4a-2 の例外はそのまま通します。
- **照合の表示:** `run.json` の `comparison` は常に `not_matched`、`physical_validity` は常に `not_established` です。検証で確かめたのは、与えた核と模型の定義に対する次数ごとの畳み込みの数値です。模型の物理の妥当性は確かめていません。実測スペクトルとは比べていません。
- **範囲の外:** 実測スペクトルとの比較（M4b-3）、斜出射、弾性散乱、軌道ごとの長さや自己無撞着な長さ、窓の外の損失、内因性の損失の分離、装置関数、背景の差引き、フィットは扱いません。QUEELS-XPS の完全な実装ではありません。

[数式と単位（E8）](docs/EQUATIONS.md)・[M4b の設計](docs/M4B_DESIGN.md)・[M4b-1 の事前登録](docs/M4B1_PREREG.md)・[検証記録](docs/VALIDATION.md)

### 束縛エネルギーの入出力（開発版 0.4.0.dev0、M4b-2）

**公開した 0.3.0 には入っていません。** `multiple_loss_spectrum` のキーワード引数 `energy_axis` と、`xps-multiple-loss` の `--spectrum` の新しいスキーマ `xpsloss-measured-spectrum-v0` で、一次線と出力の端をフェルミ準位基準の束縛エネルギーの軸で与え、結果を同じ軸で受け取り、測定の条件を記録できます。物理は足しません。計算は M4b-1 の模型のまま運動エネルギーの軸で行い、束縛エネルギーは入口と出口の変換 E = C − B だけです。`energy_axis` を与えなければ（既定）、結果は M4b-1 とビット単位で同じです。`multiple_loss_order_areas` と単一損失の関数は `energy_axis` を取りません。`examples/measured_spectrum.json` は合成の入力です。スキーマの名前に measured とありますが、実測のデータではなく、測定の記録の値もすべて例のために作ったものです。

```powershell
.\.venv\Scripts\python.exe -m xpsloss xps-multiple-loss examples/synthetic_bulk.json --spectrum examples/measured_spectrum.json --loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total --model first_loss_13a_chain_v0 --max-order 3 --epsrel 1e-7 --out runs/multiple-loss-binding
```

```python
import numpy as np
from xpsloss import (BindingDeltaPrimary, BindingEnergyAxis, EnergyAxisReference, Material, Oscillator,
                     multiple_loss_spectrum)

material = Material("synthetic_single_oscillator_not_real_material", 0.0, (Oscillator(15.0, 180.0, 2.0, 0.5),))
axis = BindingEnergyAxis(  # synthetic values, not a real measurement
    photon_energy_eV=1253.5,
    reference=EnergyAxisReference(kind="photon_and_offset", kinetic_energy_offset_eV=4.25),
    primary_reference_binding_eV=249.25,
    primary_reference_choice="synthetic delta position, not a measured peak",
    source={"name": "synthetic source, not a real instrument", "monochromatic": True,
            "satellite_treatment": "none (synthetic)"},
    calibration="not_recorded",
    charge_correction={"applied": False, "method": None, "shift_eV": None},
    core_level="synthetic level",
    emission_angle_deg=0.0,
    acceptance_half_angle_deg="not_recorded",
)
primary = BindingDeltaPrimary([249.25], [1.0])
edges = np.arange(248.25, 309.75, 1.0)

result = multiple_loss_spectrum(material, 1000.0, primary, edges, loss_window_eV=(1.0, 50.0),
                                weight_length_A=20.0, epsilon_momentum="k_total",
                                model="first_loss_13a_chain_v0", max_order=3, epsrel=1e-7, energy_axis=axis)
info = result.energy_axis
print(info["conversion_constant_eV"], info["kinetic_output_edges_eV"][[0, -1]], info["not_recorded"])
print(info["support_order_binding"], info["support_order_kinetic"], result.truncation["state"])
print(result.spectrum[:3], sorted(result.outside))
```

2つの例は、M4b-1 の例と同じ合成の入力を、C = 1249.25 eV で束縛エネルギーの軸に移したものです。`epsrel=1e-7` も同じです。249.25 eV の δ の一次線は 1000 eV に写ります。この線は内部のビンの境の上にあり、低い束縛エネルギーの側のビンに入ります。

- **換算:** `reference.kind` で2つの形のどちらかを選びます。`photon_and_offset` は、`photon_energy_eV` と `kinetic_energy_offset_eV` から C = hν − φ とします。`paired_anchor` は、同じ校正の目盛りの上の対応点 `binding_anchor_eV`・`kinetic_anchor_eV`（線の名前と出典つき）から C = B⋆ + E⋆ とします。2つの形は排他です。`paired_anchor` の含意されるずれ hν − C は記録するだけで、計算に使いません。各座標は、入力の数から C − B までを厳密に計算し、最後に1回だけ最も近い倍精度の数に丸めます。丸めの最大は `rounding_max_eV` に残ります。
- **光子エネルギーは必須です。** `photon_energy_eV` は、どちらの形でも数値で与えます。線源の名前から数値を補いません。`"not_recorded"` にもできません。
- **測定の記録:** `source`、`calibration`、`charge_correction`、`core_level`、`emission_angle_deg`、`acceptance_half_angle_deg` は必須の欄ですが、`"not_recorded"` を書けます（対応点の線の名前と出典、`primary_reference_choice` も同じ）。そのときは警告 `energy_axis_field_not_recorded` を出し、欄の名前を `energy_axis["not_recorded"]` に並べます。記録は結果と `run.json` に写すだけで、計算に使いません。観測の計数を入れる欄はありません。
- **パッケージは軸をずらしません。** 入力は、すでに校正・補正の済んだフェルミ準位基準の目盛りにあるとみなします（警告 `binding_axis_assumed_calibrated` を常に出します）。記録した校正のずれ、帯電の補正の量、出射角は「何をしたか」の記録です。仕事関数や帯電の量を推測しません。真空準位を基準にした束縛エネルギーは扱いません。`electron_energy_eV`（E0）はこれまでどおり明示し、束縛エネルギーから決めません。C − B_ref との差は記録に残ります。
- **拒否:** 負の束縛エネルギーと、減少の順の端は受け付けます。運動エネルギーが負になる座標は `invalid_energy_axis` で拒否します。一次線の位置・一次線のビンの端・出力の端の全部を合わせた異なる座標が、丸めの後に順序を保たなければ、`axis_conversion_not_injective` で拒否します。どちらも `MultipleLossInputError` で、核を評価する前に出ます。
- **出力の向き:** 配列は、入力の束縛エネルギーの端の向き（増加でも減少でも）で返します。`bin[i]` は常に入力の端の i 番目と i+1 番目の間のビンです。外の領域の名前は、入力の向きによらず `outside["low_binding"]`・`outside["high_binding"]` です。`multiple_loss_spectrum.csv` の端の列は `binding_edge_low_eV`、`binding_edge_high_eV`、`kinetic_edge_low_eV`、`kinetic_edge_high_eV` の4つです。内部の境の上の δ は、低い束縛エネルギーの側のビンに入ります。
- **`single_loss` は運動エネルギーのままです。** 結果に入っている M4a の結果 `result.single_loss` は変えません（運動エネルギーの昇順、`outside["low"]`・`outside["high"]`）。
- **台の判定は両側で行います。** 表示のビンに届く次数の上限を、元の束縛エネルギーの厳密な差と、丸めた運動エネルギーの差の両方で求め（`support_order_binding`、`support_order_kinetic`）、弱い側の主張を返します。`truncation.support_order` は大きい方、ビンが `complete` になるのは両側とも届かないときだけ、上界には大きい方の幅を使います。ビンの上界は、そのビンと、それより低い束縛エネルギーの側の全部に落ちる量の累積の上界です。2つの側が食い違っても拒否せず、警告 `axis_conversion_support_disagreement` を出し、件数と先頭 8 件までの例を `energy_axis["support_comparison"]` に残します。返す数値は、丸めた運動エネルギーの座標の問題の値で、換算の丸めは返す誤差に入れません。
- **照合の表示:** `comparison` は常に `not_matched`、`physical_validity` は常に `not_established` です。検証で確かめたのは、換算、向き、記録と、束縛エネルギーの実行が、独立に換算した運動エネルギーの入力での同じ実行とビット単位で同じことです。模型の物理の妥当性は確かめていません。
- **範囲の外:** 実測スペクトルとの比較（M4b-3）、斜出射、弾性散乱、軌道ごとの長さ、内因性の損失の分離、装置関数、観測の計数の入力は扱いません。QUEELS-XPS の完全な実装ではありません。

[数式と単位（E8a）](docs/EQUATIONS.md)・[M4b の設計](docs/M4B_DESIGN.md)・[M4b-2 の事前登録](docs/M4B2_PREREG.md)・[検証記録](docs/VALIDATION.md)

文献由来の固定入力と実行例は [benchmarks/README.md](benchmarks/README.md)、出典・帰属表示・和則の検算は [docs/BENCHMARK_INPUTS.md](docs/BENCHMARK_INPUTS.md)、複素誘電応答の設計と規約は [docs/DIELECTRIC_DESIGN.md](docs/DIELECTRIC_DESIGN.md) にあります。

仮想環境のactivateは不要です。通常のインストールには依存パッケージ取得のため通信が必要ですが、計算処理そのものは外部APIを使いません。依存関係を事前に導入すればオフラインで計算できます。

## macOS / Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m xpsloss capabilities --json
.venv/bin/python -m xpsloss bulk examples/synthetic_bulk.json --out runs/bulk-demo
.venv/bin/python tests/run_suite.py
```

**`examples/` の例題は合成材料です。Al、CoO、Siなどの実測・文献データとして利用しないでください。** 文献由来の入力は `benchmarks/inputs/` にだけ置き、出典ごとに分けています。既存の出力ディレクトリには上書きしません。再計算するときは新しい出力名を指定します。

計算結果の `run.json` には、コマンドラインで指定した表記のままの入力パス `input_file` と、入力ファイルのバイト列の SHA-256 `input_sha256` が入ります。

API の数値未収束では `NumericalConvergenceError` またはその派生例外が発生します。CLI は数値未収束を `not_converged`、入力の誤りを `error` として標準エラーに返します。どちらも終了コード2で止まり、新しい出力先は作りません。`dielectric` と `xps-fixed-depth` は、成功時の標準出力に `warning_count` を含めます。`run.json` の診断名 `static_sign`、`kinematics` の旧名は[変更履歴](CHANGELOG.md)にあります。標準出力と標準エラーの JSON は、ASCII 以外の文字をエスケープして書くので、端末・パイプ・ファイルの文字コードによらず正しい JSON になります。`run.json` は UTF-8 で、文字をそのまま保存します。

## ライセンスと引用

Keisuke Sato がコードを MIT、自作文書を CC BY 4.0 で提供します。文書中の第三者の引用は各権利者の条件に従います。第三者データ（`benchmarks/sources/**` と `benchmarks/inputs/*elfdb*`）は両方の適用範囲から除外され、N. Pauly、F. Yubero、S. Tougaard による *ELF dielectric functions for various materials*, v1.0（2022）の CC BY 4.0 が適用されます。[MITライセンス](LICENSE)・[帰属表示](NOTICES.md)・[ELFデータ集](https://doi.org/10.5281/zenodo.6024064)・[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

ソフトウェアの引用には Keisuke Sato、xpsloss と、使った版を記載してください。DOI [10.5281/zenodo.24490457](https://doi.org/10.5281/zenodo.24490457) はすべての版をまとめたもので、最新版を指します。0.3.0 の DOI は [10.5281/zenodo.24575549](https://doi.org/10.5281/zenodo.24575549)、0.2.0 の DOI は [10.5281/zenodo.24543166](https://doi.org/10.5281/zenodo.24543166)、0.1.0 の DOI は [10.5281/zenodo.24490458](https://doi.org/10.5281/zenodo.24490458) です。データ集を使用した場合はその出典も引用してください。[引用情報](CITATION.cff)

## 検証と資料

- `docs/VALIDATION.md`: この配布物で実行した検証と未検証事項。
- `docs/EQUATIONS.md`: 数式、単位、検証の根拠。
- `docs/M1_EQUATION_MAP.md`: 1997年原著の式(13)–(17)の転記、記号と単位、項ごとの実装との対応、未確定の点。
- `docs/M2_FIRST_CASE.md`: M2の最初のケース（Al、1000 eV、法線出射）、参照ソフトとの比較の許容幅（比較の前に固定）。
- `docs/REVIEW_CHECKLISTS.md`: 実装前・実装後・参照比較のレビュー項目。
- `docs/QUEELS_DISPERSION_TEST.md`: 分散のある入力で、QUEELS-XPS v2.1 の出力が公表式と候補式のどちらに一致するかを照合した事前登録試験。内部の実装は未確認。
- `docs/ROADMAP.md`: 公開までのロードマップと、初版についての決定。
- `docs/NEXT_TASK.md`: 次の工程。
- `docs/M2_2_DESIGN.md`: ギャップのある材料の固定深さ（M2-2）の設計と決定。
- `docs/M2_2_PREREG.md`: M2-2 の検証の条件、許容幅、期待する結果（参照値と実装の前に固定）。
- `docs/M3A_DESIGN.md`: 法線出射の深さ平均（M3a）の設計と決定。
- `docs/M3A_PREREG.md`: M3a の検証の条件、許容幅、期待する結果（参照値と実装の前に固定）と、参照値・経路の値の凍結の記録。
- `docs/M4A_DESIGN.md`: 単一損失のスペクトルと損失窓の面積（M4a）の設計と決定。
- `docs/M4A_PREREG.md`: M4a の検証の条件、許容幅、期待する結果（参照値と実装の前に固定）と、参照値の凍結の記録。
- `docs/M4B_DESIGN.md`: 法線出射の多重損失（M4b）の設計と決定。
- `docs/M4B1_PREREG.md`: M4b-1 の検証の条件、許容幅、期待する結果（参照値と実装の前に固定）と、参照値の凍結の記録。
- `docs/M4B2_PREREG.md`: M4b-2（束縛エネルギーの入出力と測定の記録）の検証の条件、許容幅、期待する結果（参照値と実装の前に固定）と、参照値の凍結の記録。
- `docs/DIELECTRIC_DESIGN.md`: 複素誘電応答の方式比較、符号・幅の規約、ギャップ端と極の扱い。
- `docs/BENCHMARK_INPUTS.md`: 固定入力の出典、利用条件、光学和則の監査。
- `docs/PROJECT_PLAN.md`: 完全モデルとQuanty接続へ進むための設計。
- `docs/REFERENCES.md`: 一次資料と用途。
- `NOTICES.md`: ライセンスの適用範囲、同梱データの帰属表示、参照ソフトとの関係。

試験は、ソースのフォルダーで `tests/run_suite.py` を実行します。ランナーは、常にこのフォルダーの `src` を試験します。モジュールの探索先と、子プロセスに渡す `PYTHONPATH` の先頭に `src` を置き、試験するパスを表示し、読み込んだ xpsloss が `src` のものでなければ止まります。テストが0件のときも失敗します。`python -m unittest discover -s tests -v` を使うときは、`PYTHONPATH` に `src` を設定するか（PowerShell では `$env:PYTHONPATH = 'src'`、bash では `PYTHONPATH=src`）、編集可能な導入（`pip install -e .`）にしてください。そうでなければ導入済みの別の版を試験することになり、`tests/test_suite_source.py` が失敗します。試験ファイルを1つだけ指定して `unittest` で実行すると、この確認は働きません。

開発版（0.4.0.dev0、M4b-2、コミット `7a11252`）では、Windows の作業ツリーから433件中417件が合格し、16件は skip です。M4b-2 の事前登録の82行は、すべて登録した状態になり、台帳の25,948キーはすべて合格しました。このうち1,771キーは、実装担当とは別の担当が、パッケージのコードを使わずに有理数で作った参照との照合です。`energy_axis` を与えない呼び出しは M4b-1 のコードとビット単位で同じで、事前登録した14の変異はすべて検出されました。独立レビューの指摘は blocker 0、major 0、minor 5 で、1件を直し（コミット `7a11252`。上の数はこのコードでのもの）、4件は既知の点として記録しました。これらが確かめるのは換算、向き、記録で、模型の物理の妥当性ではありません。実測スペクトルとは比べていません。CI の結果は push の後に記録します。M4b-1（コミット `ac1b656`）では、Windows の作業ツリーから400件中384件が合格し、16件は skip でした。M4b-1 の事前登録の109行は、すべて登録した状態になりました。実装担当とは別の担当が作った解析の参照の3,572キー（パッケージのコードを使わない）と、ELF データ集の Al・Si の参照の6,096キー（公開の深さ平均の点の値だけを共有する別の経路）は、すべて事前登録の許容幅で合格しました。1次の欄は M4a・M4a-2 とビット単位で同じで、独立レビュー（blocker 0、major 0、minor 10）の指摘を直した前と後の2回とも、事前登録した29の変異がすべて検出されました。これらが確かめるのは畳み込みの数値で、模型の物理の妥当性ではありません。実測スペクトルとは比べていません。コミット `e7fa6a3`（CI が Linux と macOS で見つけた DS0 の全軸の残りの区間の丸めを直した版）では、CI（run 46072925415）の全5ジョブに合格しました。テストの各ジョブでは401件中385件が合格して16件が skip、Git の履歴がない展開した sdist では382件が合格して19件が skip でした。0.3.0（M4a、段階 R4 のランナーの修正、M4a-2 を含む）では、Windows の Python 3.12.10 で、作業ツリーから347件中332件が合格し、15件は skip です（求めたときだけ動く全件の検査と、行わなかった V9 の検査）。M4a-2 の事前登録の465行は、すべて登録した状態になりました。実装担当とは別の担当が、パッケージのコードを使わずに、凍結した M4a の参照の表から写像した4,176キーはすべて合格し、重みの長さと乗数の既定の選び方は M4a とビット単位で同じで、事前登録した26の変異はすべて検出されました。この写像したキーが確かめるのは、乗数・重みの長さ・分母の配線と別経路の積分で、核の値の独立な照合ではありません。M4a-2 のコミット `c20ac79`（0.3.0.dev0、同じコード）では、CI（run 45902228392）の全5ジョブに合格しました。テストの各ジョブでは347件中332件が合格して15件が skip、Git の履歴がない展開した sdist では329件が合格して18件が skip でした。段階 R4 では337件中322件が合格していました。段階 R4 のコミット `e4e024f` では、CI（run 45761092337）の全5ジョブに合格しました。テストの各ジョブでは337件中322件が合格して15件が skip、Git の履歴がない展開した sdist では319件が合格して18件が skip でした。段階 R4 の前の M4a のコミット `25bf472` では、333件中318件が合格し（15件は skip）、CI（run 45748442743）の Linux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブにも合格し、2つの照合も合格しました。0.2.0 では、Windows の Python 3.12.10 で、作業ツリーから275件中262件が合格し、13件は skip です（求めたときだけ動く全件の検査と、行わなかった V9 の検査）。コミット `30286b7`（0.2.0.dev0、同じコード）では、CI（run 45578321049）のLinux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブにも合格し、2つの照合も合格しました。M3a の前の全208件は、Windows の Python 3.12.10 で作業ツリーから、Python 3.10.19（依存の最低版）で書き出したコピーから合格しました。コミット `fd31dcd` では、CI（run 45481433276）の Linux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブにも合格し、2つの照合も合格しました。0.1.0 の全157件は、Windows の Python 3.12.10 と 3.10.22（依存の最低版）で、作業ツリーから合格しました。最後の修正（段階4f）の前の公開候補では、全119件が、Windows の Python 3.12.10（NumPy 2.5.5、SciPy 1.18.1）と Python 3.10.22（NumPy 1.24.0、SciPy 1.10.0）で、作業ツリーから合格しました。Python 3.12.10 では、展開した sdist からも合格しました。UTF-8 モードを使わず、リダイレクトした出力は既定の文字コード（cp932）でした。その前の版の全112件は、両方の環境で展開した sdist からも合格しました。別担当も、Windows の Python 3.11.14（NumPy 2.4.6、SciPy 1.17.1）で、展開した sdist の全112件を確認しました。以前の版のテスト（75件）は、CI で Linux（Python 3.10 と依存の最低版、3.12）、Windows（3.12）、macOS（3.13）と、wheel の導入・展開した sdist のテストのジョブに合格しました。コミット `9b29e01` の公開候補（157件）は、CI（run 45325087478）で Linux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブに合格しました。2つの照合も合格しました。その前の公開候補 `5625455`（131件）と `c7fdd22`（119件）も、同じジョブに合格しています。詳細は `docs/VALIDATION.md` を参照してください。
