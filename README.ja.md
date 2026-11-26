# xpsloss

公開理論に基づくXPSの誘電損失計算を、Python APIとCLIで再現するためのパッケージです。0.2.0 の対象は M0、M1、M2-1 と、硬いギャップのある材料の固定深さ（M2-2）、外から与える一定の長さでの法線出射の深さ平均（M3a）です（v0.1.0 は M0、M1、M2-1）。開発版（0.3.0.dev0、未公開）には、法線出射の単一損失スペクトルと有限の損失窓の面積（M4a）を加えています。[English README](README.md)

**現時点で動くのは、次の6つです（6つ目は開発版 0.3.0.dev0 だけ）。**
- **振動子で表した損失関数**
- **その因果的な複素誘電応答**
- **均質・等方な無限媒質のバルク非弾性散乱**
- **法線出射・固定深さでの表面・正孔・真空中の項（Simonsen 1997 の公表式。ギャップのない材料と、0.2.0 からは硬いギャップのある材料）**
- **法線出射の深さ平均（原著の式(18)、外から与える一定の長さ）**
- **法線出射の単一損失スペクトルと有限の損失窓の面積（原著の式(19)、外から与える一定の長さ。開発版）**

**斜出射、軌道ごとの長さや自己無撞着な長さでの深さ平均、多重散乱を含むスペクトルは未実装です。QUEELS-XPSの完全再現版ではありません。分散のあるAlでは、QUEELS-XPS v2.1 との比較で、両方の運動量指定が事前の基準を満たしませんでした。**

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
| 法線出射の単一損失スペクトルと有限の損失窓の面積、外から与える一定の長さ（開発版 0.3.0.dev0、M4a） | 実装済み・テストあり・別担当がレビュー。原著の式(19)のゼロ損失線と単一損失の衛星を、損失窓の中だけで計算する。事前登録の 100 行が登録した状態になり、実装担当とは別の2つの担当が作った参照の 2041 キー（解析の参照 927、ELF データ集の Al・Si の通し経路 1048、動く端の合成入力 66）と、事前登録の許容幅で一致した。成功は推定誤差が予算の内という意味で、証明された上界ではない。実測スペクトルとは比べていない |
| 軌道ごとの長さ・自己無撞着な長さでの深さ平均 | 未実装 |
| 多重散乱、窓の外の損失、内因性の損失の分離、装置関数、背景の差引き | 未実装 |
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

次の2つは開発版（0.3.0.dev0）だけにあります。`examples/single_loss_spectrum.json` は合成の一次線（1000 eV の δ 線1本）と出力の端です。

```powershell
.\.venv\Scripts\python.exe -m xpsloss xps-loss-window-area examples/synthetic_bulk.json --loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total --out runs/loss-window-area-demo
.\.venv\Scripts\python.exe -m xpsloss xps-single-loss examples/synthetic_bulk.json --spectrum examples/single_loss_spectrum.json --loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total --out runs/single-loss-demo
```

深さ平均（M3a、0.2.0 から）は、`effective_cross_section_depth_averaged` と `xps-depth-averaged` で、法線出射の固定深さの11出力を、原著の式(18)の単一散乱の重み x exp(−x/L)/L² で深さ方向に平均します。重みの長さ L（`weight_length_A`、`--weight-length-A`）は必須で、呼び出す側が与え、全経路で一定です。出所は `--weight-length-source` に書けます。深さの積分は運動量の被積分関数の中で解析的に閉じるので、深さの格子はありません。11出力は固体内の経路長当たりの Å⁻¹ eV⁻¹ で、真空中の損失を含み、`total = no_hole + hole_charge * hole` が成り立ちます。原著の留保（式(18)は REELS の重みの流用で XPS での妥当性は検討を要すること、結果は最初の1回の非弾性散乱だけを表し intrinsic の部分を含むこと）は `warnings` に書きます。軌道ごとの長さや自己無撞着な長さ、弾性散乱、一様でない深さ分布、損失窓の面積、スペクトルは扱いません。未収束は `DepthAverageConvergenceError`（`reason` つき）です。事前登録の検査で使ったのは、ELF データ集の Al・Si と合成入力、ほぼすべて 1000 eV、長さ 5、20、80 Å と、極限の確認の 10⁻⁶–10⁶ Å、両方の運動量の指定です。`run.json` は照合した条件を記しますが、実行ごとの範囲の判定はしません。[M3a の設計](docs/M3A_DESIGN.md)・[事前登録](docs/M3A_PREREG.md)

単一損失のスペクトルと損失窓の面積（M4a、開発版 0.3.0.dev0）は、`loss_window_area`、`single_loss_spectrum` と、`xps-loss-window-area`、`xps-single-loss` で計算します。原著の式(19)を、法線出射、比例定数 1、外から与える一定の長さ L で使い、ゼロ損失線 J0 = F と単一損失の衛星 J1 = L∫F(E+T)K_sc(E0,T;L)dT を返します。核は深さ平均（M3a）の出力で、すべての一次線のエネルギーに対して設定の `electron_energy_eV`（E0）に固定します。L は式(18)の重みと式(19)の乗数の両方に使います。

- **損失窓は必須です。** `--loss-window-eV a b`（`loss_window_eV`、0 < a < b < E0）に既定値はありません。原著は式(19)の積分が発散するとして下限 0.05 eV から積分しており（p. 1618）、結果は下限に依存しえます。窓の外の損失は計算せず、収支にも入れません。下限への感度（`cutoff_sensitivity`）を一緒に返します。
- **面積:** `window_area.csv` に、11列（8項、`no_hole`、`hole`、`total`）の窓の面積（Å⁻¹）、推定誤差、予算、下限への感度（無次元）が入ります。`run.json` には `single_loss_weight`（L×`total` の面積、無次元。確率ではなく 1 を超えうる）と、面積の逆数とその区間（Å。平均自由行程とは呼びません）が入ります。
- **一次線の形式:** `--spectrum` の JSON に、`schema`（`xpsloss-single-loss-spectrum-v0`）、`axis`（`kinetic_energy_eV`）、`primary`、`output_edges_eV` を書きます。`primary` は、δ 線（`kind: "delta"`、`energies_eV`、`intensities`）か、ビン（`kind: "binned"`、`edges_eV`、ビンごとの積分強度 `intensities`）のどちらか一方です。強度は非負で、位置は非負の運動エネルギーです。出力の端は配列か `{"min", "max", "step"}` で与えます。一次線は内因性の損失の衛星を含まないものとして扱います。
- **ビン積分と単位:** `single_loss_spectrum.csv` は、出力のビンごと（左閉右開、最後のビンだけ右端を含む）と両側の外（`low`、`high`）について、ゼロ損失、衛星の3成分（`total`、`no_hole`、`hole`）、その和 `spectrum` のビン積分と推定誤差を、一次線の強度の単位で返します。負の値は切りません。密度（強度の単位/eV）は `run.json` の `*_per_eV` にあります。設定ファイルの `loss_grid_eV` は検証しますが使いません。
- **成功の意味:** 成功は、面積と各ビンの推定誤差が予算の内という意味で、証明された誤差の上界ではありません。解像の根拠の強さは `diagnostics.resolution_basis` に記録します。未収束は `SingleLossConvergenceError`、入力の誤りは `SingleLossInputError`（どちらも `reason` つき）で、CLI は終了コード2で止まり、出力先を作りません。
- **照合の表示:** `run.json` の `comparison` の `matched` は、入力の条件が、配布する完了した比較の表の登録した条件（ELF データ集の Al・Si、登録した一次線・出力の端・窓・長さ・精度）と数値として完全に一致し、その必須キーがすべて合格だったことだけを意味します。核の条件だけが合えば `kernel_conditions_checked`、ほかは `not_matched` です。上の合成の例は登録した条件の外なので、`matched` にならず、警告が出ます。`physical_validity` は常に `not_established` です。
- **範囲の外:** 多重散乱、斜出射、軌道ごとの長さや自己無撞着な長さ、弾性散乱、窓の外の損失、内因性の損失の分離、装置関数、背景の差引き、フィットは扱いません。活動線のある材料で、ギャップ端をまたぐ窓とギャップ下の窓は拒否します。QUEELS-XPS の完全な実装ではありません。

[数式と単位](docs/EQUATIONS.md)・[M4a の設計](docs/M4A_DESIGN.md)・[事前登録](docs/M4A_PREREG.md)・[検証記録](docs/VALIDATION.md)

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

ソフトウェアの引用には Keisuke Sato、xpsloss と、使った版を記載してください。DOI [10.5281/zenodo.24490457](https://doi.org/10.5281/zenodo.24490457) はすべての版をまとめたもので、最新版を指します。0.2.0 の DOI は [10.5281/zenodo.24543166](https://doi.org/10.5281/zenodo.24543166)、0.1.0 の DOI は [10.5281/zenodo.24490458](https://doi.org/10.5281/zenodo.24490458) です。データ集を使用した場合はその出典も引用してください。[引用情報](CITATION.cff)

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
- `docs/M4A_DESIGN.md`: 単一損失のスペクトルと損失窓の面積（M4a、開発版）の設計と決定。
- `docs/M4A_PREREG.md`: M4a の検証の条件、許容幅、期待する結果（参照値と実装の前に固定）と、参照値の凍結の記録。
- `docs/DIELECTRIC_DESIGN.md`: 複素誘電応答の方式比較、符号・幅の規約、ギャップ端と極の扱い。
- `docs/BENCHMARK_INPUTS.md`: 固定入力の出典、利用条件、光学和則の監査。
- `docs/PROJECT_PLAN.md`: 完全モデルとQuanty接続へ進むための設計。
- `docs/REFERENCES.md`: 一次資料と用途。
- `NOTICES.md`: ライセンスの適用範囲、同梱データの帰属表示、参照ソフトとの関係。

開発版（0.3.0.dev0、M4a を含む）の全333件は、Windows の Python 3.12.10 で作業ツリーから合格しています（15件は skip。求めたときだけ動く全件の検査と、行わなかった V9 の検査）。コミット `25bf472` では、CI（run 45748442743）の Linux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブにも合格し、2つの照合も合格しました。0.2.0 の全275件は、Windows の Python 3.12.10 で作業ツリーから合格しています（13件は skip。求めたときだけ動く全件の検査と、行わなかった V9 の検査）。コミット `30286b7`（0.2.0.dev0、同じコード）では、CI（run 45578321049）のLinux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブにも合格し、2つの照合も合格しました。M3a の前の全208件は、Windows の Python 3.12.10 で作業ツリーから、Python 3.10.19（依存の最低版）で書き出したコピーから合格しました。コミット `fd31dcd` では、CI（run 45481433276）の Linux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブにも合格し、2つの照合も合格しました。0.1.0 の全157件は、Windows の Python 3.12.10 と 3.10.22（依存の最低版）で、作業ツリーから合格しました。最後の修正（段階4f）の前の公開候補では、全119件が、Windows の Python 3.12.10（NumPy 2.5.5、SciPy 1.18.1）と Python 3.10.22（NumPy 1.24.0、SciPy 1.10.0）で、作業ツリーから合格しました。Python 3.12.10 では、展開した sdist からも合格しました。UTF-8 モードを使わず、リダイレクトした出力は既定の文字コード（cp932）でした。その前の版の全112件は、両方の環境で展開した sdist からも合格しました。別担当も、Windows の Python 3.11.14（NumPy 2.4.6、SciPy 1.17.1）で、展開した sdist の全112件を確認しました。以前の版のテスト（75件）は、CI で Linux（Python 3.10 と依存の最低版、3.12）、Windows（3.12）、macOS（3.13）と、wheel の導入・展開した sdist のテストのジョブに合格しました。コミット `9b29e01` の公開候補（157件）は、CI（run 45325087478）で Linux（Python 3.10 と依存の最低版、3.12）、macOS（3.13）、Windows（3.12）と配布物のジョブに合格しました。2つの照合も合格しました。その前の公開候補 `5625455`（131件）と `c7fdd22`（119件）も、同じジョブに合格しています。詳細は `docs/VALIDATION.md` を参照してください。
