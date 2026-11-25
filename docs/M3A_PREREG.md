# M3a の事前登録　検証の条件・許容幅・期待する結果

作成日: 2026-11-25。状態: **固定（実装と参照値の計算の前）。** 根拠は [M3a の設計](M3A_DESIGN.md)（利用者が §10 の U1–U10 を推奨どおり決定、2026-11-25）の §7–§9。この文書と期待する状態の表を含むコミットの後に、参照値の計算と実装を始める。比較の後に、許容幅と期待する結果を変えない。変える場合は、この文書の旧版を残し、別の試験として記録する。不合格なら、深さ核、項の組立て、誤差の判定、参照経路の順に切り分けて記録する。

凍結の前に、Codex が案を3回確かめた（1回目: 修正後に凍結、blocker 3件、major 7件、minor 2件。2回目: 修正後に凍結、新しい blocker 2件、major 5件。3回目: 修正後に凍結、前回の blocker は解消、明文の不整合4か所）。指摘はすべて反映した（§7）。

記号は設計の「表記と状態」に従う。\(K_{{\rm sc},j}\) は式(18)の深さ平均（Å⁻¹ eV⁻¹）、\(\lambda\) は重みの長さ（Å、引数 `weight_length_A`）、\(j\) は11出力で、順序は `term_13a`、`term_13b`、`term_14`、`term_15_T1`、`term_15_T2`、`term_15_T3`、`term_16_e`、`term_16_h`、`no_hole`、`hole`、`total` である。

## 1. 入力

数値はすべて十進で与え、API には binary64 に変換して渡す。参照の計算にも同じ binary64 の値を渡す。機械で読める同じ内容が、`benchmarks/numerics/m3a/expected_states_m3a.json` の `materials` にある。分散の値は無次元、強度は eV²、他は eV。

| 記号 | 内容 | 来歴 |
| --- | --- | --- |
| Al | `benchmarks/inputs/al_elfdb_v1.0.bulk.json`（SHA-256 `978931c1651881e6c7dfa739f4d626e4496225bf24bceb815466ff651934eae2`）の材料。記録 `al_elfdb_v1.0.json`（`17ded0465a5e42124eebf44d9573369a5cbea87df59b5df11216ea2f7432d1c4`）、元データ `ELF.txt`（`1a59cc61a562f53e3bda0eec0c29df3d0a84b601fdc0aecee454be62db738f48`）。\(g=0\)、\((E,A,\Gamma,\alpha)=(15.6,243.3598,1.5,1.0)\) | ELF データ集 v1.0（CC BY 4.0） |
| Al0 | Al の分散だけを 0 にしたもの | synthetic |
| Y0、Y5 | \(g=0\)、\((18,200,3,\alpha)\)、\(\alpha=0\)、0.5 | synthetic |
| Z | Al の強度を 0 にしたもの | synthetic |
| W4、W6、W7、W8 | \(g=0\)、\((15,180,\Gamma,0)\)、\(\Gamma=3\times10^{-4}\)、\(3\times10^{-6}\)、\(3\times10^{-7}\)、\(3\times10^{-8}\) | synthetic（段階4f の線） |
| W7a | \(g=0\)、\((15,180,3\times10^{-7},2\times10^{-8})\) | synthetic |
| W2 | \(g=0\)、\((15,81,4\times10^{-7},0)\) と \((18,180,2\times10^{-7},0)\) | synthetic（`tests/test_numerics_4f.py` と同じ） |
| Si | `benchmarks/inputs/si_elfdb_v1.0.bulk.json`（`471e513416369eec3dd12302a8e00823113935f78988293205b34f469f7200be`）の材料。記録 `si_elfdb_v1.0.json`（`49020931fcb6d83f799b6a8b0cfd0a6eedb8c3a5dc446b88793bc9e4e5704079`）、元データ `ELF_Si.txt`（`109330f657a62a50ca43e1077fe74b20810cb5477b1dd1122b202de84c17de3d`）。\(g=1.12\)、\((10,6.08458,5,0.5)\)、\((14,30.4229,5,0.5)\)、\((16.8,212.9603,3.8,0.5)\) | ELF データ集 v1.0（CC BY 4.0） |
| S0 | Si の分散だけを 0 にしたもの | synthetic |
| Z0 | Si の強度をすべて 0 にしたもの | synthetic |
| Z1 | Si に強度 0 の線 \((0.5,0,0.3,0.5)\) を加えたもの | synthetic |
| N1 | \(g=2\)、\((3,6,0.02,0.5)\) | synthetic（M2-2 と同じ） |
| N2 | \(g=2\)、\((3,6,0.2,0.5)\)。N1 の線幅を10倍にした、ギャップより上の成功対照 | synthetic（設計 §8.2、U10） |
| P1、P0、M1、UE、SP | [M2-2 の事前登録](M2_2_PREREG.md) §1 と同じ | synthetic |

## 2. 損失の配列、精度の組、予算

**損失の配列**（SHA-256 は binary64 のリトルエンディアンのバイト列）:

| 記号 | 式（NumPy） | 点数 | SHA-256 |
| --- | --- | --- | --- |
| GA | `0.2*np.arange(1, 500)` | 499 | `6214a4e35ef20fa9f1c9f5f19a376d6888847f081e5942672eba0ae9b9abf33e` |
| GC | `np.linspace(0.2, 99.8, 499)`（CLI が設定から作る格子） | 499 | `e69297cdd297d6432c4976f071938e792fb93f29f8681cba47da75677a842524` |
| F | `0.2/4*np.arange(4, 1997)`（F[::4] は GA とビット単位で等しい） | 1993 | `0c204e56a600e9d92793f98b5747494da583c0c841528f7ba65f7fc03765fe47` |
| F2 | `F[::2]` | 997 | `a493a70b1103fd493ce43ea0dd56a2734224955a9f4cc791a2fec6309484221f` |
| FRD | `np.r_[F[::-4], F[:40:4]]`（逆順と重複） | 509 | `9a33b993c2304b0bf4a4a3f0201f361b9f1b17a99972660973f7ddb772fcbe3c` |
| GA10G | `np.r_[GA[:10], 1.12]` | 11 | `b220dc4f8bd483ebe3751b9b736eb14b98e3b9ec510160524d9714234b866909` |

GA と GC は 499 点中 340 点が binary64 で異なる。単点の損失は、表の JSON に十進と binary64 の16進表記の両方で書いた。ギャップ下の Si と S0 の単点には GA の最初の5点を使う（GA の3点目は `float("0.6")` と異なる）。

**精度の組** \((\varepsilon_{\rm abs}\) [Å⁻¹ eV⁻¹]、\(\varepsilon_{\rm rel}\)、limit\()\):

| 記号 | 値 | 渡し方 |
| --- | --- | --- |
| DEF | \((10^{-12},10^{-9},200)\) | API は明示。CLI は設定の `integration` に書く |
| TIGHT | \((10^{-14},10^{-11},200)\) | 同上 |
| IMP | \((10^{-30},10^{-15},200)\) | 同上 |
| LOOSE | \((10^{-6},10^{-9},200)\) | 同上 |
| LIM1 | \((10^{-12},10^{-9},1)\) | 同上 |
| APIDEFAULT | DEF と同じ値 | API に `epsabs`、`epsrel`、`limit` を渡さない（関数の既定値の検査、U6） |
| CLIDEFAULT | \((10^{-9},10^{-7},200)\) | CLI の設定に `integration` を書かない（設定の既定値の検査、U6） |
| CLIPREC | CLIDEFAULT と同じ値 | API に明示する（CLIDEFAULT と同じ精度の API の対照） |

**共通の条件:** 電子のエネルギーは 1000 eV（MS2 は 4000 eV、UE は 3 eV）、出射角 0°（RJ4 の1件だけ 10°）。

**予算:**

- 実装が自分の返値で判定する予算: \(B^{\rm run}_j=\max(\varepsilon_{\rm abs},\varepsilon_{\rm rel}|\widehat K_j|)\)。
- 外部の照合で参照値から作る予算: \(B^{\rm ref}_j=\max(\varepsilon_{\rm abs},\varepsilon_{\rm rel}|K_{j,\rm ref}|)\)。精度の組は、その実行のもの。参照の区間 \(|K_{j,\rm ref}|\pm u_j\) から作る下側と上側を \(B^{\rm ref-}_j\)、\(B^{\rm ref+}_j\) とする。
- \(e_j\) は、実装が公表する出力 \(j\) の推定誤差（`component_errors[j]`）。

## 3. 期待する状態

**表:** `benchmarks/numerics/m3a/expected_states_m3a.json`（SHA-256 `d87dad8f02a2fd211679229bd7c5205d908dbb2b3b241b30a27e94469f414ed5`、UTF-8、LF）。生成は `python benchmarks/numerics/m3a/expected_states_m3a.py`（NumPy と標準ライブラリだけを使い、`xpsloss` を読み込まない）。`--check` で再生成してバイト単位で比べる（それ以外の引数は拒否する）。

**数:** 2041 の実行、66965 の点のキー。区分ごとに MS 66377、NC 43、EI 20、RJ 525。参照のキーは 498、V4・V5 の極限のキーは 192。

**照合の予定:** 表の `check_plan` に、検査ごとの（実行の ID、元の添字、参照のキー）の組を列挙した。参照値は物理のキーで共有し、合否は「検査、実行の ID、元の添字、出力の列」の単位で数える。

| 検査 | 内容 | 全件の組 | 定例の組 |
| --- | --- | --- | --- |
| V7 | 参照のキーを持つ実行の各点（成功した EI を含む） | 1582 | 1020 |
| V3(d) | MS8、MS9、MS14、MS15 の 13a と 20 Å の参照 | 384 | 288 |
| V4 | MS8、MS14 と極限のキー | 264 | 168 |
| V5 | MS9、MS15 と極限のキー | 120 | 120 |
| V11(iii) | 活動線のある材料の、ギャップ下の単位電荷の MS と EI の全点（参照のキーが一致する点だけ (1) を行う） | 249 | 101 |
| V11(v) | ギャップ下の単位でない電荷の MS の全点（参照、Z0 は解析的な 0） | 424 | 129 |

実行の前に予定の組の数を確かめ、実行の後に、合格、不合格、照合未確定、EI の指定の未収束、未実行の内訳を数え、予定の数と照合する。欠けた組を合格に数えない。

**生成のときの検査**（`assert` ではなく例外で止める。`python -O` でも働く）: 実行のキーの重複がない。点のキーが2つの区分を持たない。拒否以外の実行の点がすべてその実行の区分である。拒否の添字が存在し、範囲内で、重複しない。EI の実行は1点で、許す理由が §4 の語彙の空でない部分集合である。NC の理由が語彙にある。実行が指す参照のキーと 13a の照合のキーがすべて定義されている。EI の実行はすべて参照のキーを持つ。さらに、族の定義とは別に書いた表と照らして、EI の族ごとの許す理由、添字を必須とする拒否（RJ1、RJ2、RJ2-SPedge、RJ3）と添字のない拒否（RJ2-SPsign、RJ4、RJ5）、参照のキーの有無の物理の規則（MS5、MS6、MS11、MS13、MS16、MS17、EI は全点、MS14 はギャップ下、MS15 は 1 eV）、ギャップ下の単位でない電荷の点に参照があること、V4・V5 の予定が極限のキーを持つこと、実行が指す各参照のキーの材料・エネルギー・指定・長さ・損失が、その実行のその元の添字のものと一致すること（13a の照合のキーは長さ 20 Å）を確かめる。生成の入力をわざと誤らせた4例（EI2 に直接和の理由を足す、RJ2 の添字を消す、MS14 の参照を外す、Y5 の 5 Å と 80 Å の参照先を入れ替える）は、いずれもここで止まった。

**キー:** 点のキーは、材料、電子のエネルギー、運動量の指定、正孔電荷（binary64）、重みの長さ（binary64、または不正な値の種類）、損失（binary64）、精度の組、出射角、入口（API か CLI か）、CLI の出所の文字列、出力先の事前の状態と名前（ASCII か日本語か）である。実行のキーは、さらに損失の順序つきの列を含む。拒否の実行では `offending_indices` の点だけをキーに数える。添字を持たない拒否（RJ4 の不正な長さと角度、RJ2-SPsign、RJ5）は、入力の全点を数え、評価は0回とする。

**CLI の実行:** 各実行の `cli_config` に、設定ファイルの `loss_grid_eV` と `integration`（精度の組の名前、または書かない `null`）を書いた。試験は、入力の設定ファイル（材料は表のとおり）の `loss_grid_eV` と `integration` をそれに置き換え、実行の前に設定を読み戻して、損失のビット列と精度が表と一致することを確かめる。

| 実行 | `loss_grid_eV` | `integration` |
| --- | --- | --- |
| MS4（CLI）、RJ5 | `{"min": 0.2, "max": 99.8, "step": 0.2}`（GC） | DEF（MS4 の1件は書かない、CLIDEFAULT） |
| NC2 | `{"min": 5.0, "max": 5.0, "step": 1.0}`（単点 5 eV） | IMP |
| RJ3 | `{"min": 1.999, "max": 1.999, "step": 1.0}`（単点 1.999 eV） | DEF |

出力先の名前は、`output_dir_name` が `ja` の実行で `深さ平均の出力`（表の `cli_strings`）、他は ASCII。RJ5 の既存の出力先には、内容が既知のファイルを1つ置いておき、実行の後にそのバイト列の SHA-256 が変わらないことを確かめる。

**区分:**

| 区分 | 期待 | 判定 |
| --- | --- | --- |
| MS 必須成功 | 例外なしに返り、全点で11出力の推定が有限・非負で \(e_j\le B^{\rm run}_j\) を満たす | 表の `expect` の追加の条件（厳密な 0、直接和の診断、警告の有無）も確かめる |
| RJ 特定理由の拒否 | 本体の計算の前に、指定の型と理由で止まる | 型、`reason`（その属性を持つ型）、違反の添字、本体のパネルの評価が0回 |
| NC 特定理由の未収束 | `DepthAverageConvergenceError` で止まる | `reason`、`loss_eV`、`index`、指定があれば `column` と `subdivisions` |
| EI 参照に合う成功、または指定の未収束 | どちらかになる | 成功なら MS と同じ判定と、§6 の V7 の3値の照合。失敗なら `reason` が許す集合にあり、`index` が 0 |

CLI の実行は、状態（`complete`、`not_converged`、`error`）、終了コード、出力先の事後の状態（作られた、作られない、内容が変わらない）で判定する。stderr の JSON は型ごとに確かめる: `not_converged` は `error_type`、`reason`、損失と添字。`FixedDepthGapDomainError` は `error_type`、`reason`、違反の損失と添字。`FileExistsError` は `error_type` と文言だけ（損失、添字、理由を持たない）。

`expect` の `direct_sum_diagnostics` が `null` のときは、`diagnostics["direct_sum"]` が長さ \(n\) の全 `null` のリストで、`direct_sum_max` が `null` である。`nonnull` のときは、少なくとも1点が数の辞書である。

**族の要約**（各実行の具体的なキーは表にある。根拠は設計 §8.2）:

| 族 | 区分 | 内容 |
| --- | --- | --- |
| MS1 | MS | Al、Al0、Y0、Y5 × 両指定 × \(Q_h\in\{0,0.5,1\}\) × \(\lambda\in\{5,20,80\}\) Å × GA |
| MS2 | MS | Al、4000 eV、両指定、\(Q_h=1\)、20 Å、GA |
| MS3 | MS | Al、Si、`k_total`、\(Q_h=1\)、20 Å × F、F2、FRD（V10(ii)） |
| MS4 | MS | `k_total`、\(Q_h=1\)、20 Å、GC。API: Al、Si、Z（DEF）と Al（APIDEFAULT、CLIPREC）。CLI（新しい出力先）: Al（出所なし、DEF）、Al（日本語の出所と日本語の出力先、DEF）、Al（CLIDEFAULT）、Si（直接和の診断が空でない）、Z（低損失の警告なし） |
| MS5 | MS | Al、Al0、Y0、Y5 × \(T\in\{0.05,0.2,1,5,11,15.6,30,60,99.8\}\) の単点 × 両指定 × \(Q_h\) 3値 × \(\lambda\) 3値。V7 の参照点 |
| MS6 | MS | W4、W6、W7 × \(T=\sqrt{135}\)（binary64）× 両指定 × \(Q_h\in\{0,1\}\) × 20 Å。V7 の参照点 |
| MS7 | MS | Z × GA × 両指定 × \(Q_h\in\{0,1\}\) × 20 Å と、Z × 5 eV × `k_total` × IMP。11出力が厳密に 0 |
| MS8、MS9 | MS | Al、Y0 × \(T\in\{1,15.6,30\}\) × 両指定 × \(Q_h\) 3値 × \(\lambda\in\{10^{-6},10^{-5}\}\)（V4）、\(\{10^5,10^6\}\)（V5） |
| MS10 | MS | Si、S0 × GA × 両指定 × \(Q_h\) 3値 × \(\lambda\) 3値 |
| MS11 | MS | Si、S0 × {GA の最初の5点、\(g\pm10^{-3}\)、\(g\pm10^{-6}\)、1.2、2、5、11、16.8、17、30、60、99.8} の単点 × 両指定 × \(Q_h\) 3値 × \(\lambda\) 3値。V7 の参照点 |
| MS12 | MS | Z0 × GA10G × 両指定 × \(Q_h\in\{0,1\}\) × 20 Å（端を含めて厳密に 0、ギャップ由来の警告なし）。Z1 × GA × 両指定 × \(Q_h=1\) × 20 Å（Si の MS10 と同じ点が受理され、ギャップ下の各点の証明が成り立つ） |
| MS13 | MS | N2 × 4 eV × 両指定 × \(Q_h\in\{0,1\}\) × 20 Å。直接和の診断が `null`。V7 の参照点 |
| MS14 | MS | Si、S0 × GA の最初の5点 × 両指定 × \(Q_h\in\{0,0.5\}\) × \(\lambda\in\{10^{-6},10^{-5}\}\)。同じ点で \(Q_h=1\) と LOOSE。\(T\in\{1.2,5,16.8\}\) × \(Q_h\) 3値 × 同じ \(\lambda\)。V4。ギャップ下の点は V7・V11(v) の参照点でもある |
| MS15 | MS | Si × \(T\in\{1,1.2,5,16.8\}\) × 両指定 × \(Q_h\) 3値 × \(\lambda\in\{10^5,10^6\}\)。V5。1 eV は V7・V11(v) の参照点でもある |
| MS16 | MS | Al × MS5 の \(T\)、Si × MS11 の \(T\) の単点 × 両指定 × \(Q_h\in\{-0.5,1.5\}\) × 20 Å。V6 |
| MS17 | MS | Si × GA の最初の5点 × 両指定 × \(Q_h=\) `nextafter(1,2)` × 20 Å × LOOSE。直接和の診断が `null`（V12(b)）。V7 の参照点 |
| NC1 | NC | Al、Si × 5 eV × `k_total` × \(Q_h=1\) × 20 Å × IMP。`panel_floor_exceeds_budget`、`subdivisions` 0 |
| NC2 | NC | NC1 の Al を CLI で。`not_converged`、終了コード2、出力先を作らない、`reason` が同じ |
| NC3 | NC | Si、S0 × GA の最初の5点 × 両指定 × \(Q_h=1\) × \(\lambda\in\{10^{-6},10^{-5}\}\)。`direct_sum_floor_exceeds_budget`、列 `total`、添字 0、`subdivisions` 0 |
| EI1 | EI | Si × `k_total` × \(Q_h=1\) × 1 Å の単点 0.2、0.4 eV。`budget_not_met`、`panel_floor_exceeds_budget`、`direct_sum_floor_exceeds_budget` |
| EI2 | EI | Al × `k_total` × \(Q_h=1\) × 20 Å × TIGHT の単点 0.05、0.2 eV。`budget_not_met`、`panel_floor_exceeds_budget` |
| EI3 | EI | W8 × \(\sqrt{135}\) × 両指定 × \(Q_h\in\{0,1\}\)（DEF）。W7 × \(\sqrt{135}\mp1.5\times10^{-7}\)、`k_parallel`、TIGHT。W2 × 12 eV × `k_parallel` × TIGHT。W7a × \(\sqrt{135}\) × 両指定（DEF）。いずれも 20 Å。`budget_not_met`、`panel_floor_exceeds_budget`、`surface_denominator_unresolved` |
| EI4 | EI | N1 × 4 eV × 両指定 × \(Q_h\in\{0,1\}\) × 20 Å。`budget_not_met`、`panel_floor_exceeds_budget`、`surface_resonance_unresolved`、`d3_ill_conditioned`（M2-2 の固定深さでは、この点は D3 の悪条件で止まった） |
| EI5 | EI | Al × `k_total` × \(Q_h=1\) × 20 Å の単点 0.01、0.02、0.1 eV。`budget_not_met`、`panel_floor_exceeds_budget`。各点の状態を V16 に記録する |
| RJ1 | RJ | Si × GA10G × 両指定（`DielectricGapEdgeError`、`offending_indices` [10]）。Si × \([g+\delta_g/2]\)（[0]） |
| RJ2 | RJ | P1 × 1.999 eV と 1.9 eV、P0 × 1.999 eV（`surface_pole_not_excluded`）、M1 × 2.5 eV（`monotonicity_condition_failed`）、いずれも両指定。UE × 3 eV × `k_parallel` × 1.5 eV（`surface_pole_not_excluded`、`details["points"][0]` に `U_eV` の欄があり `None`）。`FixedDepthGapDomainError` の `details["points"][*]["index"]` が [0]、`offending_loss_eV` がその損失 |
| RJ2-SPedge | RJ | SP × 1.0 eV × 両指定。`DielectricGapEdgeError`、[0]（端の検査が静的符号より先に止める） |
| RJ2-SPsign | RJ | SP × 2.5 eV × 両指定。静的符号の `ValueError`（文言に `static sign tolerance`。損失の添字なし） |
| RJ3 | RJ | P1 × 1.999 eV × `k_total` を CLI で（単点の設定）。`error`、終了コード2、出力先を作らない。stderr の `error_type` が `FixedDepthGapDomainError`、`reason` が `surface_pole_not_excluded`、違反の添字 [0] |
| RJ4 | RJ | Al × 5 eV × `k_total`。長さ `0.0`、`-1.0`、`nan`、`inf`、`1+0j`、`True`（この型で渡す）は `ValueError`（文言に `weight_length_A`）。\(\theta=10^\circ\) は `ValueError`（文言に `theta_deg`）。LIM1 は `ValueError`（文言に `integration limit too small` と `loss_eV=5, index=(0,)`、既存の固定深さと同じ規則） |
| RJ5 | RJ | MS4 の Al の CLI と同じで、出力先が既存（日本語の名前）。`error`（`FileExistsError`）、終了コード2、既知のファイルの内容が不変 |

**定例の試験と全件の記録:**

- 表の各実行に `routine` を付けた（1329 実行）。定例の試験（CI を含む）は、この部分集合を実行する。Si の順序・重複の配列（MS3 の FRD）、NC3 の両方の長さ（GA[0]）、Si の 5 Å の低損失（MS11、GA[0]、\(Q_h=1\)）を定例に含めた。
- 全件の記録は、表の全 2041 実行と、§6 の各検査の全件である。統括担当が1回実行し、結果、所要時間、扱える精度の下限の観察を [検証記録](VALIDATION.md) に書く。
- 実行時間の都合で定例の部分集合を狭めるときは、区分と許容幅と全件の対象を変えずに、差分を記録する。
- 事前検査の受理や NC3 の床の下限だけから、全件の成功を結論しない。

| 検査 | 定例の試験 | 全件の記録 |
| --- | --- | --- |
| MS の予算、RJ・NC・EI の状態（V13） | `routine` の実行 | 全実行。区分ごとに件数を集計する |
| V1、V2 | 全点 | 同じ |
| V3(d)、V4、V5、V7、V11(iii)・(v) | `check_plan` の組のうち `routine` の実行のもの（上の表） | `check_plan` の全組。合格、不合格、照合未確定を別々に数える |
| V3(a)–(c)、V6 | `routine` の実行の参照点 | 全参照点 |
| V8 | 必須の3点（凍結した経路の値と尺度） | §6 の全点 |
| V9 | \(\lambda=20\) Å | 全点 |
| V10 | (i) は \(\lambda=20\) Å、(ii) は Al と Si の FRD（比べる相手は同じ条件の MS1・MS10 の GA） | (i) は全長さ、(ii) は MS3 の全配列 |
| V11(i)・(ii)・(iv)、V14 | 全件 | 同じ |
| V15 | MS4 の CLI と API の対照、RJ5、`run.json` の読み戻しは全件。11列の診断の再計算は定例の MS の 1255 実行 | 同じ対照と、全 MS の 1951 実行の再計算 |
| V12 | — | レビュー担当が §6 の変異の表の全件を1回行う |
| V16 | — | 記録 |

## 4. 実装が守る約束（試験が使う面）

**API**（新しいモジュール `xpsloss.xps_depth_averaged`。関数、結果の型、例外を `xpsloss` から公開する）:

```python
effective_cross_section_depth_averaged(
    material, electron_energy_eV, loss_eV, *, weight_length_A, epsilon_momentum,
    theta_deg=0.0, hole_charge=1.0, epsabs=1e-12, epsrel=1e-9, limit=200,
) -> DepthAveragedResult
```

- `DepthAveragedResult` の欄: `electron_energy_eV`、`weight_length_A`、`theta_deg`、`hole_charge`、`epsilon_momentum`、`loss_eV`、8項、`no_hole`、`hole`、`total`（いずれも長さが損失の数の1次元配列）、`quadrature_error`（各損失の11推定の最大）、`component_errors`（11出力の名前から1次元配列への辞書。ギャップの有無によらず必須）、`diagnostics`（辞書）、`warnings`（文字列のタプル）。`depth_A` と `path_length_A` は持たない。
- `diagnostics` の必須の欄:
  - `direct_sum`: 長さ \(n\) のリスト。ギャップ下の直接積算を使った点 \(i\) は辞書 `{"K_dir", "e_dir", "delta_rec", "e_total"}`、使わなかった点は `None`。辞書の4つの数は有限で、誤差は非負。`e_total` は `component_errors["total"][i]` と等しく、`delta_rec` は binary64 で計算した `abs(total[i] - K_dir)` と等しい。
  - `direct_sum_max`: 使った点だけの `e_dir` と `delta_rec` の最大（辞書）。該当する点がなければ `None`。
  - `error_budget_ratios`: 11出力の名前から、\(e_j/B^{\rm run}_j\) の1次元配列への辞書。
  - `max_budget_ratio`: `{"value", "column", "index"}`。同率なら、出力の順序が先の列、次に元の添字の小さい点を選ぶ。全推定が 0 なら `{"value": 0.0, "column": "term_13a", "index": 0}`。
  - ギャップのある材料では、固定深さの `FixedDepthResult.diagnostics` と同じ名前の欄（ギャップの事前検査の結果、ギャップ下の各点の証明、`max_subgap_identity_residual`）。
- **深さ核**（V1 が使う）: `xpsloss._depth_kernel.depth_averaged_weights(k, b, p, *, e_k=0.0, e_b=0.0, e_p=0.0)`。\(k\)、\(b\)、\(p=1/\lambda\) は Å⁻¹ の binary64 のスカラー、`e_*` はその入力の絶対誤差。返値は辞書で、キー `W13`、`N`、`W14`、`W15b`、`W15c`、`W16e`、`W16h` のそれぞれが、尺度化した重み \(p^2W/(2b)\)（設計 §3.1 の \(PW/M\)）の `(re, im, err_re, err_im)`。`err_*` は、入力の誤差と演算の丸めを含む上界である。符号つきの分子の相殺（零交差）では、値と誤差を返す（例外にしない）。尺度化した値そのものが binary64 の指数の範囲を外れ（正規数としても非正規数としても表せず）0 か無限大になるときだけ、`FloatingPointError` を送出する。本体は同じ核を使い、この例外を `weight_out_of_range` にする。
- **本体のパネルの評価:** 深さ平均の運動量の被積分関数の評価は、初回のバッチ、個別の再試行、細分のいずれも `xpsloss.xps_depth_averaged._panel_integrals` を通す。試験はこれを差し替えて呼出しを数える。拒否では0回である。成功する対照の実行で、差し替えた関数が呼ばれることも確かめる。
- **V10(i) の試験用の口:** `xpsloss.xps_depth_averaged._depth_averaged_point(material, electron_energy_eV, loss_eV, *, weight_length_A, epsilon_momentum, hole_charge, epsabs, epsrel, limit, panels=None, adapt=True)` が、1つの損失の11出力の値、推定、最終のパネルの端（\(\ln k\) の昇順の列）を返す。`panels` を与えて `adapt=False` なら、そのパネルで再適応せずに評価する。
- **固定深さの関数**（`effective_cross_section_fixed_depth`）の引数、返値、失敗は変えない（V14）。

**失敗:**

- 数値の未収束は `DepthAverageConvergenceError(NumericalConvergenceError)`。属性: `reason`（下の表の文字列）、`loss_eV`（float）、`index`（元の配列の 0 始まりの int）、`column`（出力の名前か `None`）、`estimate`、`budget`（有限の float か `None`。11予算の判定に達した後の失敗でだけ数）、`subdivisions`（止まるまでに行った細分の回数、int）、`panels`（止まった時点のその損失のパネルの数、int）、`inner`（内側の原因の文言か `None`。既存の固定深さの例外に由来するときは、その文言を含む）。
- **段階の順**（M2-2 の順序を引き継ぐ）: (1) 入力の検査 → (2) 格子全体の端の検査 → (3) 運動量の区間 → (4) 格子全体の単調性と静的符号 → (5) ギャップ下の各点の証明 → (6) 線の位置と初回のパネル → (7) 初回のバッチの節点の評価 → (8) 初回のバッチの後の床の判定（直接和の床を、正のパネルの床より先に判定する） → (9) 細分。拒否は既存の型と属性を保つ。端の拒否（`offending_indices`）と `FixedDepthGapDomainError`（`details["points"][*]["index"]` と `offending_loss_eV`）は、既存の場所に違反したすべての点の添字を返す（例えば P1 の [1.999, 1.9] は [0, 1]）。静的符号と入力の `ValueError` は、登録した型と文言で判定し、添字を求めない。(3) で区間が決まらない入力は、(4) より前に `momentum_bounds_unresolved` で止まる。`DepthAverageConvergenceError` は1点だけを報告し、最も早い段階で止まった点のうち元の添字の小さい点を選ぶ。

| `reason` | 使う条件 |
| --- | --- |
| `panel_floor_exceeds_budget` | 初回のバッチの後で、細分でも係数の精細化でも減らない誤差だけで、ある出力の \(B^{\rm run}_j\) を超える。対象は、(i) 初回のバッチで決まる正のパネルの丸めの床と、(ii) 同じパネルで係数を精細化（設計 §3.1 の選択的な再計算）した後にも残る係数と重みの構成の誤差の積分、の和である。精細化で減りうる初回の係数の誤差だけでは、この理由にしない。細分しない。`estimate` はその床、`budget` は予算（NC1） |
| `direct_sum_floor_exceeds_budget` | ギャップ下の直接和の床 \(f_{\rm dir}\)（設計 §3.1 の下限 \({\rm eps}\cdot\int\Sigma_{\rm dir}\,d\ln k\) 以上）が `total` の予算を超える。初回のバッチの後に判定し、細分しない（NC3） |
| `budget_not_met` | パネルの数が `limit` に達するまで細分しても、ある出力の推定が予算を超える（固定深さの「did not converge with N panels」に当たる）。`panels` は `limit`、`estimate` \(>\) `budget` |
| `surface_denominator_unresolved` | 節点で表面の分母が解像できない（固定深さの「surface denominator unresolved at quadrature node」） |
| `oscillator_denominator_unresolved` | 節点で振動子の分母が解像できない（「oscillator denominator unresolved at quadrature node」） |
| `linewidth_unresolved` | 線幅や線の位置が運動量のパネルで解像できない（「linewidth unresolved on momentum panels」、「surface oscillator sum unresolved」、「unresolved oscillator and surface lines」） |
| `surface_resonance_unresolved` | ギャップありの表面の共鳴の探索が決着しない（`_gap_response` の surface search の失敗） |
| `d3_ill_conditioned` | D3 の悪条件（`_D3ConditionError` に由来するもの） |
| `coefficient_unresolved` | D3 の悪条件以外の係数の失敗（残存帯、除去帯の導関数の区間、減衰の解像） |
| `momentum_bounds_unresolved` | 運動量の区間が決まらない |
| `weight_out_of_range` | 深さ核が `FloatingPointError` を送出した（尺度化した重みが指数の範囲を外れる。零交差の相殺は含まない） |
| `response_out_of_range` | 誘電応答や振動子の和が数値の範囲を超える |

- **失敗の側の観測**（NC と、EI が失敗したとき）: 試験は `_panel_integrals` を包み、呼出しごとに、その損失の元の添字、評価したパネルの組（\(\ln k\) の左右の端）、返った各パネルの11出力の値と推定を記録する（`_panel_integrals` は、1回の呼出しで1つの損失のパネルの組を受け、形（パネルの数, 11）の値と推定を返す。11出力の値は、同じパネルの組の値の和に、集約の丸めだけを加えたものとする）。例外の値を、この記録と照合する。
  - 床の理由: `index` と `loss_eV` が記録の損失と一致する。`column` が違反の列。`budget` が、記録した最後のパネルの組の値の和から作った \(\max(\varepsilon_{\rm abs},\varepsilon_{\rm rel}|\sum|)\) と一致する（相対 \(4u_r\) 以内）。`estimate` \(>\) `budget`。`subdivisions` が 0 で、記録したパネルの組が初回のものだけ（同じパネルでの係数の再評価は許す）。
  - `budget_not_met`: `panels` が記録した最後のパネルの組の数と一致し、`limit` に等しい。`estimate` が、その列の記録したパネルの推定の和以上で、`budget` を超える。
  - 節点の理由（分母、D3）: 例外が `_panel_integrals` の呼出しの中で送出された。表面の探索の理由: 例外が `xpsloss._gap_response.resolve_surface` の呼出しの中で送出された（試験はこれも包む）。`inner` が対応する既存の文言を含む。
  - 理由の名前と添字だけでは合格にしない。
- **実効の精度:** `diagnostics["effective_settings"]` に、計算に使った `{"epsabs", "epsrel", "limit"}` を書く（API と CLI の両方）。
- 入力の不正と定義域の拒否: `ValueError`（文言は §3 の RJ4）、`FixedDepthGapDomainError`（M2-2 と同じ属性）、`xpsloss.dielectric.DielectricGapEdgeError`（`offending_indices` を持つ）、静的符号の `ValueError`。CLI ではどれも `status=error`、終了コード2。
- CLI の未収束は、stderr の JSON に `status=not_converged`、`error_type`、`reason`、損失と添字、終了コード2。どの失敗でも、計算が済むまで出力先を作らない。既存の出力先は内容まで保つ。

**CLI:**

```powershell
xpsloss xps-depth-averaged CONFIG --weight-length-A 20 [--weight-length-source TEXT] --epsilon-momentum k_total [--theta-deg 0] [--hole-charge 1] --out DIR
```

- 損失の配列と精度は設定ファイルから取る（`loss_grid_eV` から `np.linspace`、`integration`、省略時は設定の既定値）。
- 出力: `xps_depth_averaged.csv`（列 `loss_eV,total,no_hole,hole`、8項、`quadrature_error`。見出しは ASCII）と `run.json`（UTF-8）。stdout と stderr の JSON は ASCII のエスケープ。
- `run.json` の欄は設計 §9 のとおり。`model_id` は `simonsen1997_depth_average_external_v0`。`weight` は `{"form": "x*exp(-x/L)/L^2", "normalization": 1, "length_A", "length_kind": "external_constant", "length_source"（文字列か null）, "depth_integration": "closed_form_laplace", "elastic_scattering": false, "depth_distribution": "uniform_semi_infinite", "loss_window_for_length": null}`。`diagnostics` に、API と同じ名前で `component_errors`、`error_budget_ratios`、`max_budget_ratio`、`direct_sum`、`direct_sum_max` を書く（配列はリスト、`None` は null）。実効の精度、入力のバイト列の SHA-256、列の単位、`warnings`、`applicability` を含む。

## 5. 独立な参照値と独立経路（実装担当とは別の担当）

担当は3つに分ける。実装（Codex）、多倍長の参照（担当 R1）、独立経路と期待する状態の照合（担当 R2）。R1 と R2 は実装のコードと結果を見ずに作り、自分の式の自己照合だけで本番との一致を宣言しない。

**R1 多倍長の参照**（V1、V3–V7、V11(iii)(1)・(v)）:

- パッケージ `xpsloss` の KK、閉形式、求積、係数、運動学の関数を使わない。ELF（E1）、硬いギャップの KK（E4a の主値）、式(13)–(17)の各項、深さの閉形式を、原論文の転記（[式の対応表](M1_EQUATION_MAP.md) §4.1）から自分で書く。定数は文書の値（\(\hbar^2/2m=3.80998212\) eV Å²、\(a_0=0.529177210903\) Å）を使う。mpmath で作業精度 50 桁以上。
- \(K_{\rm sc}\) を2通りに計算する。(a) 自分で導いた深さの閉形式。(b) 自分の**固定深さの**被積分関数（深さの因子を平均せずに持つもの）を、位相で区切って深さ方向に求積したもの。(b) では (a) の平均済みの重みを使わない。電子群と正孔群は、8項の和とは別に原式から組み立てる。
- **参照点:** 表の `reference_keys` の 498 キー（材料、エネルギー、指定、長さ、損失）。各キーで8項（単位電荷の成分）、`no_hole`、`hole`、電子群、正孔群と、それぞれの不確かさ。証明済みのギャップ下では、\(Q_h=1\) の `total` に解析的な 0 と \(u=0\) を使い、生の8項から作った残差を別の欄に保存する。電荷 \(Q_h\) の `total` は \(K_e+Q_hK_h\)、不確かさは \(u_e+|Q_h|u_h\) とする（直接に評価してもよい）。**保存:** 値は十進 30 桁以上の文字列で保存し、各値の \(u\) に、計算値から保存した十進への丸めの幅を含める。試験は保存した十進を有理数として厳密に読み、電荷倍と群の和も有理数で作る（binary64 の \(Q_h\) も正確な有理数として使う）。
- **極限のキー:** 表の `limit_keys` の 192 キー。`lambda_to_0` では \(G_j(0)\)（設計 §2.2 の1次元求積）、包絡 \(\int(|sc_{sj}|+2k|c_{dj}|)\,d\ln k\)、\(\kappa_j\)。`lambda_to_inf` では \(G_{\infty,j}\)、登録した各長さでの過渡の包絡 \(E_{\rm tr}\)、一次の係数 \(A_1\) と一次の後の残りの上界。いずれもそれぞれの求積の不確かさ付き。
- **V1:** §6 の V1 の点の入力（binary64 の16進）。各点で、7つの重みの実部・虚部と \(n\) を、binary64 の入力の正確な有理数から**厳密な有理数**で計算し、比較用の値とする（保存は十進 40 桁以上の文字列と、有理数の分子・分母）。別に、8つの固定の点で、固定深さの因子を位相で区切って深さ方向に求積した値（変数 \(t=pa\)、相殺を避ける式か十分な追加の桁、尾と次数を変えた不確かさを記録）を作り、厳密な値との一致を導出の確認として記録する。入力の箱の8隅の厳密な値。4つの変異（\(p+2k\to p+k\)、\(s\leftrightarrow\bar s\)、`weight_14` の係数 2、符号）を自分の式に入れ、少なくとも1点で V1 の幅を超えて離れることを記録する。零交差の直後の点で、公表の誤差が \(e_n\) の寄与だけで正になる（他の寄与は \(|\hat n|=0\) に比例して消える）ことを、有理数で確かめて記録する。
- **計算の手順:** ギャップより上の主値の特異点の扱い、細い共鳴の両側の分割、長さと電荷で共有する中間量（KK と係数は材料・エネルギー・指定・損失ごとに一度）を、生成スクリプトに書く。代表の1条件で (a) と (b) の所要時間を測って報告する。
- **参照の不確かさ:** (a) と (b) の差、50 桁と 80 桁の差、次数・尾・位相の分割を変えた差を、キーと列ごとに保存し、合成の式を書く。そのキーを使う実行のうち最も厳しい精度の組と電荷で \(u_j\le0.01B^{\rm ref}_j\) を満たすことを確かめる（例えば EI2 の点は TIGHT）。満たさない点は細かくして計算し直す。決着しない点は「照合未確定」と列挙する。
- **凍結:** `benchmarks/numerics/m3a/generate_references_m3a.py`、`benchmarks/numerics/m3a/references_m3a.json`（入力、条件、各出力、\(u_j\) とその内訳、計算法、ライブラリの版、所要時間。ASCII、LF）、比較の試験 `tests/test_m3a_reference.py`。参照値の SHA-256 をこの文書の付録に記録してから、実装と比べる。
- **照合未確定の扱い:** 比較で未確定になったキーは、比較する担当から R1 に、キーと必要な精度だけを渡す。R1 は細かくする規則、旧値、新しい値と、新旧の SHA-256 を記録する。未確定を合格の件数に含めない。

**R2 独立経路と期待する状態の照合**（V2、V8、V9、V13 と MS の判定）:

- **経路1（数値の深さ求積）:** 公開の固定深さ関数だけを使い、設計 §3.2 の1のとおりに書く。変数 \(t=a/\lambda\)。固定深さの関数は1つの深さで損失の配列を受けるので、同じ深さの節点で \(\mathcal P\) の全損失をまとめて計算する。
  - **積分域:** \(t\in[0,r]\)、基本は \(r=40\)。区間の端は、登録した 0、0.01、0.04、0.12、0.3、0.7、1.4、2.5、4.5、8、14、25、40 に、次の2つを加えたもの（いずれも \([0,r]\) の内部だけ）。被積分関数は \(w(x)K(x)=\lambda^{-2}e^{-x/\lambda}(\kappa x+G(x))\) で、\(G\) は \(e^{-kx}\)（\(k_{\parallel-}\le k\le k_{\parallel+}\)）と \(e^{\mp ibx}\) の重ね合わせである。振動 \(e^{\mp ibx}\) は必ず \(e^{-kx}\) と組で現れるので、損失 \(T\) の振動は \(t\ge t_c(T)=40/(k_{\parallel-}(T)\lambda)\) で \(e^{-40}\) 未満に減衰する（\(k_{\parallel-}\) は式(17)の \(k_\parallel\) の下端。小さい損失では \(b\) より小さく、\(t_c\) が \(r\) を超えうる。\(b\) より大きいのは全運動量の下端である）。(i) 幾何の端 \(r2^{-m}\)（\(m=1,2,\ldots\)、\(r2^{-m}\ge1/(8k_{\parallel+,\max}\lambda)\) の範囲。\(k_{\parallel+,\max}\) は \(\mathcal P\) の損失の上端の最大）。(ii) 各区間を、その区間の中で位相の増分 \(b_{\rm eff}(t)\lambda\,\Delta t\le\pi/4\) になるまで等分する。\(b_{\rm eff}(t)\) は、\(t<t_c(T)\) を満たす \(\mathcal P\) の損失の \(b(T)\) の最大（なければ分割しない）。各区間に Gauss–Legendre 10 点。同じ深さの節点で \(\mathcal P\) の全損失を計算する。
  - **変える規則**（1つずつ変え、基本との差を列ごとに記録する）: (A) 次数 16 点。(B) 位相の増分 \(\le\pi/8\)、幾何の端の比 \(\sqrt2\) で分けた別の節点集合。(C) \(r=50\)（\([40,50]\) を同じ分割の規則で加える）。
  - **尾:** \(t>r\) の寄与の上界 \(|\kappa_j|(1+r)e^{-r}+C_je^{-r}/\lambda\)。\(C_j\) は \(|G_j(a)|\) の全深さの包絡で、設計 §2.2 のモード係数の絶対値の積分 \(\int(|c_{0j}|+|c_{sj}|+|c_{dj}|)\,d\ln k\) とその求積の誤差の上側から作る。係数の \(S\)、\(R\) には公開の誘電応答の関数を使ってよい（尾の上界にだけ使い、\(e^{-r}\) 倍で効く。値の経路ではない）。包絡を作れない点は、経路の照合を保留する。
  - \(u_j=\max(|\Delta_A|,|\Delta_B|,|\Delta_C|)\)＋節点の誤差の正の重み付きの和＋尾の上界。
  - **集合の固定:** V8 の点の集合 \(\mathcal P\)（§6 の V8）を、材料、エネルギー、指定、\(Q_h\)、\(\lambda\) ごとに先に固定する。予備の計算で \(\mathcal P\) 全体の \(\max|R_{\rm total}|\) と節点の精度の配分を決め、再評価で最終の \(D_j\)、\(u_j\) と、全節点で \(\alpha_ie_{ij}\le f_i\eta D_j\)（\(\eta=1/9\)、\(f_i\) は等分）を確かめる。定例の試験でも、この凍結した尺度を使う。
  - ギャップなしでは、節点の共通の `quadrature_error` を各列に保守的に配る。他の正孔電荷の `total` を \(Q_h=1\) の節点の成分から作るときは、誤差を \(e_e+|Q_h|e_h\) と演算の丸めで配り、ギャップ下の単位電荷の直接和の小さい誤差を流用しない。
  - 次数、打切り、位相で区切った別の節点集合、尾の包絡、節点の誤差の正の重み付きの和を、列ごとに記録する。未取得の点、必須の3点の失敗、保留をそれぞれ記録する。
- **経路3（I3 型）:** 凍結した `tests/xps_route_b_reference.py` の \(k_\perp\) の核を、経路1の基本の規則と同じ深さの節点で平均する（`k_parallel`）。変える規則は、経路1の (A)(B)(C) と、(D) `kpar_nodes` の `n_base` を 12 から 24 にしたもの、(E) I3 のモジュールの \(k_\perp\) の中心と尾の固定の分点の配列を倍に細かくしたもの（(D)(E) はモジュールの関数の引数や定数を別のプロセスで差し替えて行う。M2-2 の参照と同じ方式。凍結したファイルは変えない）。尾は経路1と同じ包絡を使う。1軸ずつ変え、各点・各列の変化の最大が V9 の幅の 1/10 以下であることを確かめる。
- **凍結:** 経路1と経路3を計算するスクリプト `benchmarks/numerics/m3a/generate_paths_m3a.py`（上の規則をそのまま書く）と、その値、尺度、不確かさとその内訳、節点の配分の `benchmarks/numerics/m3a/paths_m3a.json` を作り、JSON の SHA-256 をこの文書の付録に記録してから、実装と比べる。定例の試験は凍結した値を読む。
- **期待する状態の照合:** `tests/test_m3a_expected_states.py`（表の `routine` の実行）と、全件の記録用の `benchmarks/numerics/m3a/run_expected_states_m3a.py`（全実行）。表の JSON の SHA-256 を確かめてから使う。

## 6. 検査と許容幅

設計 §8.3 の V1–V16 を、次の数値と分担で固定する。「V7 の条件」は V7 の行の3値の判定を指す。

| 番号 | 検査 | 合格条件 | 書く担当 |
| --- | --- | --- | --- |
| V1 | 深さ核 | **点** \((k,b,p)\)（Å⁻¹）: \((10^{-3},0.5,1)\)、\((0.05,0.03,0.05)\)、\((5,0.03,0.05)\)、\((0.01,2,0.05)\)、零交差の \((0.5,1,\)`nextafter(1.5,0)`\()\)、\((0.5,1,1.5)\)、\((0.5,1,\)`nextafter(1.5,inf)`\()\)、アンダーフローの \((1,1,10^{82})\)、Al・1000 eV・\(T\in\{1,15.6,99.8\}\) eV・\(\lambda\in\{5,20,80\}\) Å の \(b\) と式(17)の \(k_\parallel\) の両端とその幾何平均（R1 が文書の定数から計算して binary64 に丸め、凍結する）。**値**（実部・虚部の成分ごと、比較の値は R1 の厳密な有理数）: 0 にならない部分は \(|\widehat W-W_{\rm exact}|\le10^{-13}|W_{\rm exact}|\)。零交差のある \(\mathrm{Re}W_{13}\)、\(\mathrm{Re}W_{15c}\) は \(10^{-13}(kL+b^2)/(pDH)\)、\(\mathrm{Im}W_{16h}\) は \(10^{-13}(kL+b^2)/(DH)\) の絶対幅（尺度化した返値と比べるので、幅に \(p^2/(2b)\) を掛ける）。\(\mathrm{Im}N\) は厳密に 0。\((1,1,10^{82})\) で尺度化した \(\mathrm{Re}W_{14}\) が 0 でなく、\(1.0000000000000002\times10^{-164}\) と相対 \(10^{-13}\) で一致。零交差の直後の点は値と誤差を返す（例外にしない）。**範囲外の対照:** \((1,1,10^{200})\) では尺度化した \(\mathrm{Re}W_{14}\approx10^{-400}\) が binary64 で表せないので、`FloatingPointError` を送出する。**誤差:** 入力の誤差 0 で、全点・全成分で \(|\widehat W-W_{\rm exact}|\le\) 公表の誤差。零交差の直後の点で、零交差のある3つの成分の公表の誤差が厳密な値の大きさ以上（この点では \(e_n\) の寄与だけが残る）。\(e_k,e_b,e_p\) をそれぞれ入力の 4 ulp にしたとき、最初の4点と零交差の直後の点の5点で、箱の8隅の厳密な値との差が公表の誤差以内（伝播の式の回帰の検査。箱の内部の上界は、実装が演算ごとの導出で保証し、報告に書く）。**変異:** R1 の記録で、4つの変異がいずれも少なくとも1点で値の幅を超える。本番の核への変異は V12 でレビュー担当が試す | R1 |
| V2 | 重みの規格化（経路1の求積器） | 経路1と同じ深さの求積器に試験関数を入れる。\(\lambda\in\{5,20,80\}\) Å、各区間 12 点。(i) \(r=8\) で再規格化の前の累積を \(1-(1+r)e^{-r}\)（定数）、\((1-e^{-r})/\lambda\)（\(1/x\)）と比べる。(ii) \(h\in\{0.01,0.1\}\)（\(t\) の値）で区間を分けて測った \([0,h]\) の \(1/x\) の寄与を \((1-e^{-h})/\lambda\) と比べる。(iii) \(r=40\) までの全積分で、\(1,\ x,\ 1/x,\ \kappa+I/x,\ e^{-sx}/x\)（\(\kappa=2\)、\(I=3\) Å、\(s=0.1+0.05i\)、\(1+0.5i\) Å⁻¹）の平均を、\(1,\ 2\lambda,\ 1/\lambda,\ \kappa+I/\lambda,\ [\lambda(1+\lambda s)]^{-1}\) と比べる（\(r=40\) の尾 \((1+r)e^{-r}\) は幅に比べ無視できる）。無次元化（定数と、すでに無次元の \(\kappa+I/x\) はそのまま、\(x\) の平均は \(\lambda\) で割り、\(1/x\) と \(e^{-sx}/x\) は \(\lambda\) を掛ける）した複素の差の絶対値が \(10^{-12}\) 以下 | R2 |
| V3 | 13a | (a) Al0、Y0、S0 の参照キーで、解析的な \(L(0,T)\log(q_+/q_-)/(\pi a_0E_0)\)（\(q_\pm\) はバルクの全運動量の端。S0 は \(T>g\)）。(b) Al、Y5、Si の `k_total` の参照キーで、R1 のバルクの別の積分。(a)(b) とも V7 の条件。(c) 証明済みのギャップ下で `term_13a` がビット単位で 0。(d) 13a は深さによらないので、MS8、MS9、MS14、MS15 の各点の `term_13a` を、表の `term_13a_reference_keys`（同じ材料・エネルギー・指定・損失の 20 Å の参照キー）の 13a と V7 の条件で比べる。分散ありの `k_parallel` にはバルクとの一致を課さない | R1 |
| V4 | \(\lambda\to0\) | MS8、MS14 のキーで、\(\widehat H_j=\lambda\widehat K_{{\rm sc},j}\) とその誤差 \(\lambda e_j\) から \(|\widehat H_j-G_j(0)|\le\lambda(|\kappa_j|+{\rm 包絡})+\lambda e_j+u_{\rm lim}+\lambda u_{\rm env}+\lambda u_\kappa\)。\(\kappa_j\) が 0 でない列（13a、`no_hole`、`total`）では R1 の \(\kappa_j\) の不確かさ \(u_\kappa\) を足し、他の列は \(\kappa_j=0\)、\(u_\kappa=0\)。13a の脱落は V3(d) で検出する。NC3 は §3 の区分で判定する | R1 |
| V5 | \(\lambda\to\infty\) | MS9、MS15 のキーで \(|\widehat H^{\rm tr}_j-G_{\infty,j}|\le E_{\rm tr}+u_{\rm tr}+\lambda\sum_{m\ne13a}|w_{jm}|e_m+\rho^H_j+u_\infty\)（eV⁻¹）。\(\rho^H_j=4u_r\,N_j\,\lambda\sum_{m\ne13a}|w_{jm}\widehat K_m|\)、\(N_j\) は和の項の数。**返却列の組立て:** \(|\widehat K_j-\sum_mw_{jm}\widehat K_m|\le\rho^K_j=4u_r\,N_j\sum_m|w_{jm}\widehat K_m|\)（13a を含む、Å⁻¹ eV⁻¹。\(\lambda\) を掛けない）。証明済みのギャップ下の単位電荷の `total` は \(|\widehat K_{\rm total}|\le\min(e_{\rm total},\varepsilon_{\rm abs})\)。**比:** \(|A_1|\) が R1 の不確かさより大きく、両方の長さで、一次の後の残りの上界と数値誤差と R1 の不確かさの**和**が \(0.1|A_1|/\lambda\) 以下の点・列にだけ、短い側と長い側の偏差の比が \([90/11,\ 110/9]\)（分数で比べる）にあることを課す。課した件数と、課さなかった件数（0 件を含む）を記録する。過渡部が厳密に 0 の列には課さない | R1 |
| V6 | ゼロ応答と線形性 | MS7 と MS12 の Z0 で11出力がビット単位で 0。MS5、MS11（20 Å、\(Q_h\) 3値）と MS16（\(-0.5\)、1.5）で、\(|\widehat K_{\rm total}-(K_{e,\rm ref}+Q_hK_{h,\rm ref})|\) を V7 の条件で（\(u=u_e+|Q_h|u_h\)）判定する。\(Q_h=0\) でも正孔の3項と `hole` を V7 で比べる | R1 |
| V7 | 多倍長の参照との照合 | 参照キーを持つ全実行（MS5、MS6、MS11、MS13、MS14 と MS15 の参照点、MS16、MS17 と、成功した EI）の各点・各出力で、\(d=|\widehat K_j-K_{j,\rm ref}|\)（有理数で厳密に計算する）、参照の不確かさ \(u_j\) について、**合格** \(d+u_j\le\min(B^{\rm ref-}_j,e_j)\)、**不合格** \(\max(0,d-u_j)>\min(B^{\rm ref+}_j,e_j)\)、その間は**照合未確定**（§5 の手順で参照を細かくする）。参照が解析的な 0 の量（ギャップ下の単位電荷の `total`、13a）は \(u=0\)、\(B^{\rm ref}=\varepsilon_{\rm abs}\)。参照がないキーは「参照未作成」と記録し、完了扱いにしない。CSV の共通の最大だけとは比べない | R1 |
| V8 | 経路1 | \(D_j=10^{-7}|R_j|+10^{-8}\max_{T\in\mathcal P}|R_{\rm total}|\)（§5 の凍結した尺度）。経路に \(u_j\le D_j/3\)、照合に \(|\Delta_j|+u_j\le D_j\)。**点の集合** \(\mathcal P\)（材料、エネルギー、指定、\(Q_h\)、\(\lambda\) ごと）: Al と Y0 の 1000 eV は MS5 の9損失、Si は MS11 の18損失、いずれも両指定、\(Q_h\) 3値、\(\lambda\) 3値。Al の 4000 eV は MS2 の GA[4]、GA[77]、GA[149]、GA[498]、両指定、\(Q_h=1\)、20 Å。**必須の点:** Al、1000 eV、`k_total`、\(Q_h=1\) の (80 Å、99.8 eV) と (20 Å、1.0 eV、15.6 eV)。経路1がこれらで成功しなければ V8 を合格にしない | R2 |
| V9 | 経路3（`k_parallel`） | Al と Al0、\(\lambda\in\{5,20\}\) Å、MS5 の9損失、\(Q_h\) 3値で、11列が \(|\Delta_j|\le10^{-6}|K_{j,\rm I3}|+10^{-6}\max|K_{\rm total,I3}|\)（最大は同じ材料・長さ・\(Q_h\) の9損失）。§5 の各規則を1軸ずつ細かくしたときの経路の変化の最大が、この幅の 1/10 以下 | R2 |
| V10 | 収束と点値の再現 | (i) `_depth_averaged_point` で、Al（MS1）と Si（MS10）の `k_total`、\(Q_h=1\)、\(\lambda\in\{5,20,80\}\)、GA の各損失。基準の実行の返値から \(B^*_j=\max(10^{-12},10^{-9}|\widehat K^{(0)}_j|)\) を決め、(a) 精度を半分にした適応の実行で各 \(e_j\le B^*_j/2\)（満たさなければさらに締め、それでも満たさなければ判定不能と記録）。(b) (a) の最終パネルを全部二分して `adapt=False` で評価し、各 \(e_j\le B^*_j/2\)。(a) と (b) の返値の差が \(B^*_j\) 以下。(ii) MS3 で、同じ binary64 の \(T\) の返値の差が、2つの公表推定の和以内 | 実装 |
| V11 | ギャップあり | (i) 事前検査の結果と理由が、同じ損失の配列での固定深さのものと同じ。(ii) §3 の RJ1、RJ2、MS12 のとおりで、拒否では本体のパネルの評価が0回。(iii) 表の `check_plan["V11iii"]`（活動線のある材料の、証明済みのギャップ下で \(Q_h=1\) の MS と EI の全点。MS3、MS4、MS10、MS11、MS12 の Z1、MS14 の LOOSE、MS15、EI1 を含む 249 組）について、(1) 参照のキーが一致する点（損失の16進表記と元の添字で結合する。例えば Si の MS3 は 39 点中 20 点、MS4 は 10 点中 8 点）で、生の8項と2群が V7 の条件で参照に合う、(2) \(|\widehat K_{\rm total}|\le\min(e_{\rm total},\varepsilon_{\rm abs})\)、(3) §4 の `direct_sum` の等式（`e_total` が `component_errors["total"]` と等しい、`delta_rec` が \(|{\tt total}-K_{\rm dir}|\) と等しい、\(e_{\rm total}=e_{\rm dir}+\delta_{\rm rec}\) が相対 \(4u_r\)）、(4) 診断の最大は該当する点だけで取り、なければ `None`。(iv) 端の3つの区別（設計 §6）を記録する。(v) 表の `check_plan["V11v"]`（証明済みのギャップ下で \(Q_h\ne1\) の MS の全点。配列で渡す MS10 の 120 点を含む 424 組）で、`total` が参照の \((1-Q_h)K_e\) と V7 の条件で合う（Z0 は解析的な 0）。各点の `direct_sum` は `None` | 実装（i–iv）、R1（iii(1)、v） |
| V12 | 変異の検出 | 下の変異の表の各行を、レビュー担当が本番のコードに1つずつ入れ、変異なしでの合格と、変異後に指定の検査が不合格になることを、同じ条件で記録する。実装は (d)–(h) の単体試験を書き、報告に試験の名前を書く | 実装、レビュー |
| V13 | 拒否と未収束 | 表の RJ、NC、EI の実行を、その `expect` のとおりに判定する（定例と全件は §3）。NC と、失敗した EI では §4 の「失敗の側の観測」も確かめる。部分的な CSV を出さない | R2 |
| V14 | 固定深さの不変 | (a) 既存の `tests/test_m22_migration.py`（v0.1.0 の `0167e38` とのギャップなしの40ケース）が変えずに合格する。(b) main の `fd31dcd` と、下の入力で、固定深さの全体が同じ環境でビット単位で一致する。**API**（同じ引数を両方の版に渡す）: (a) の40ケース、Si（10 Å、両指定、\(Q_h\in\{0,1\}\)、GA）、S0（10 Å、`k_total`、\(Q_h=1\)、GA）、Z0（10 Å、GA10G、両指定、\(Q_h=1\)）、Z1（10 Å、`k_total`、\(Q_h=1\)、GA）、N2（10 Å、4 eV、両指定、\(Q_h\in\{0,1\}\)）、N1（10 Å、1 eV、`k_total`、\(Q_h=1\)）、P1 と P0（10 Å、4 eV、`k_total`、\(Q_h=1\)）、UE（3 eV、10 Å、2.5 eV、両指定、\(Q_h=1\)。成功する）。全欄（配列はバイト列、`component_errors`、`diagnostics` と `warnings` は浮動小数を16進にした正規の JSON）。**失敗の対照**（型、文言と、その型が持つ公開の属性のすべて。例えば `reason`、`offending_loss_eV`、`details`、`offending_indices`、`gap_eV`、`exclusion_eV`、`remedy`、`diagnostics`）: N1 の 10 Å・4 eV（両指定。D3 の悪条件で、例外が `diagnostics` を持つ）、P1 の 1.999 eV、M1 の 2.5 eV、UE（3 eV、`k_parallel`、1.5 eV）、SP の 1.0 eV と 2.5 eV、Si の GA10G。**CLI**（同じ設定のバイト列を両方の版に渡す）: Si（GC、10 Å、`k_total`、DEF）、Z0（単点 min=max=1.12、10 Å、両指定、DEF）、P1（単点 1.999 eV、`k_total`）。CSV のバイト列、`run.json` と stdout・stderr の JSON（時刻、所要時間、絶対パスの欄を除く）、終了コード。**版番号**（`run.json` の `package_version`）は両方の版の比較から外し、それぞれの実行を生成したパッケージの `xpsloss.__version__` と一致することを別に確かめる（U9 の版の更新で落ちないように）。**履歴:** Git の作業ツリーでは `fd31dcd` を履歴から取り出す。CI（環境変数 `CI`）で履歴がなければ飛ばさずに失敗する。Git の作業ツリーでない展開した sdist では、理由を書いて skip する（既存の `tests/test_m22_migration.py` と同じ） | 実装 |
| V15 | CLI と記録 | MS4 で、CLI と、同じ GC と同じ実効精度の API（DEF の CLI は DEF の API、CLIDEFAULT の CLI は CLIPREC の API）の値と CSV が、ビット単位で一致する。CLI の入口から計算の関数に渡った `epsabs`、`epsrel`、`limit` を差し替えで捕まえ、設定と `run.json` の実効の値と一致することを確かめる。APIDEFAULT の実行が DEF の実行とビット単位で一致する。`run.json` を読み戻し、入力のバイト列の SHA-256、全列の値と単位、`weight`（`length_source` は日本語の実行で元の文字列、未指定で `null`）、モデル ID、実効の引数、`diagnostics` の5つの欄の型とパス（§4）、`warnings`（Z では低損失の警告なし）、`applicability`、`capabilities` の `xps_depth_averaged` を確かめる。APIDEFAULT と CLIDEFAULT の実行で、`diagnostics["effective_settings"]` がそれぞれ \((10^{-12},10^{-9},200)\)、\((10^{-9},10^{-7},200)\)。**その実行集合（定例か全件、§3）の全 MS の実行で**（R2 の照合）、`error_budget_ratios` を返値から \(e_j/B^{\rm run}_j\) として再計算して一致し、`max_budget_ratio` の値・列・元の添字が §4 の同率の規則で再計算したものと一致し、`quadrature_error` が11推定の最大と一致する。日本語の出力先の成功の実行の stdout と、RJ5（日本語の既存の出力先）の stderr で、生のバイト列が ASCII で、復号した出力先の文字列が元の日本語と一致する。RJ5 で、既知のファイルの SHA-256 が変わらず、出力先の項目が増えない | 実装、R2（全 MS の再計算） |
| V16 | 低損失の記録 | Al、`k_total`、\(Q_h=1\)、20 Å で、EI5 の3点、MS5 の 0.05 eV、MS1 の 0.2 eV から \(TK_{\rm sc}^{\rm total}\) を記録する。合否に使わない | 統括 |

**V12 の変異の表**（凍結の前に固定する。番号は設計 §8.3 の V12 の記号）:

| 番号 | 変異 | 不合格になる検査 |
| --- | --- | --- |
| K1 | 重みの \(J=p+2k\) を \(p+k\) にする | V1 の \((0.05,0.03,0.05)\) と \((5,0.03,0.05)\) の W13、W14、W15b、W16e。V7 の MS5（Al、20 Å、`k_total`、15.6 eV）の `term_13b` |
| K2 | 重みで \(s\) と \(\bar s\) を取り違える | V1 の \((0.05,0.03,0.05)\) の W13 と W16e の虚部、W16h の実部（W16h は \(-\overline{W_{16h}}\) になり、虚部は変わらない）。V7 の MS5（Al、20 Å、`k_total`、15.6 eV）の `term_13b` |
| K3 | `weight_14` の係数 2 を 1 にする | V1 の全点の W14。V7 の MS5（Al、20 Å、15.6 eV）の `term_14` |
| K4 | W15b の符号を変える | V1 の全点の W15b。V7 の MS5（Al、20 Å、15.6 eV）の `term_15_T2` |
| K5 | 前因子の \(\lambda^2\) を \(\lambda\) にする | V7 の MS5（Al、5 Å と 80 Å、15.6 eV）の 13a 以外の全列。V4 の MS8（Al、15.6 eV） |
| K6 | 16-e の定数項を落とす | V1 の全点の W16e。V4 の MS8（Al、15.6 eV）の `term_16_e`。V7 の MS5 の `term_16_e` |
| K7 | 誘電応答の運動量の指定を取り違える | V7 の MS11（Si、20 Å、5 eV、両指定）と MS5（Al、20 Å、15.6 eV、両指定）の `no_hole`。V9 |
| K8 | 13a を極端な長さで 0 にする | V3(d) の MS9（Al、\(10^6\) Å、15.6 eV）と MS15（Si、\(10^6\) Å、5 eV） |
| K9 | `total` の組立ての電荷を誤る | V6 の MS16（Al、20 Å、15.6 eV、\(Q_h=1.5\)） |
| (g) | \(e_n\) を落とす | V1 の零交差の直後の点 |
| (f) | 積の誤差の交差項 \(e_{Z_r}e_{W_r}\) を落とす | 実装の単体試験（\(\widehat Z_r=\widehat W_r=0\)、\(e_{Z_r},e_{W_r}>0\)） |
| (e) | 13b の2つの射影の一方の誤差を落とす | 実装の単体試験（実部か虚部の一方だけが誤差を担う入力） |
| (d) | パネルの集約の丸めを落とす | 実装の単体試験（複数のパネルの集約） |
| (h) | ギャップより上の係数の誤差の伝播の項（E4b、除去帯、E1 の虚部）を1つずつ落とす | 実装の単体試験（その項だけに既知の正の変化を与え、推定がその分だけ増える） |
| (a) | \(\delta_{\rm rec}\) を記録しても \(e_{\rm total}\) に足さない | 実装の単体試験（\(e_{\rm dir}<B<e_{\rm dir}+\delta_{\rm rec}\) で拒否）、V11(iii)(3) |
| (a′) | \(K_{\rm dir}\) に \(K_{\rm out}\) を写す（\(\delta_{\rm rec}=0\)） | 実装の単体試験（初回の混在のバッチと細分の両方で、丸めた8項の節点の和から \(K_{\rm dir}\) を別に作り、差が 0 に潰れない制御した入力で比べる） |
| (b) | 許されない点で直接和を使う | MS10 の \(Q_h=0,0.5\)、MS13、MS17 の診断が `None`。RJ2 の評価0回 |
| (c) | 端の上下が混在するバッチや細分で、該当しない点に直接和を使う | 実装の単体試験 |
| F1 | 直接和の床の判定を外す | NC3 の `reason` |
| F2 | 正のパネルの床の判定を外す | NC1、NC2 の `reason` |
| X1 | 深さ核の切り出しで固定深さの値を変える | V14(b) |
| X2 | 全実行で診断の数（`error_budget_ratios`、`max_budget_ratio`）をすべて 0 にする | V15 の全 MS の再計算 |
| X3 | EI と NC で、パネルを評価せずに `budget_not_met` を返す | V13 の失敗の側の観測（`panels` と評価の回数） |

単体試験で検出する行（(a)、(a′)、(c)–(h)）の制御した入力と列は、実装担当が報告に書き、レビュー担当が変異を入れる前に固定する。

## 7. 凍結の前の確認（統括担当と Codex、2026-11-25）

- **ギャップの事前検査:** 受理の条件は深さによらない（設計 §6）ので、M2-2 のコード（main `fd31dcd`）の事前検査だけを、表のギャップのある材料の実行の損失と指定で実行した（本体のパネルは評価していない。M3a の値ではない）。Si と S0 の GA と MS11 の全単点、Si の F、F2、FRD、GC、Z0 の GA10G、Z1 の GA、N2 と N1 の 4 eV、MS14・MS15・MS17・NC3・EI1 の単点は受理された。RJ1（添字 [10] と [0]）、RJ2（P1、P0 は `surface_pole_not_excluded`、M1 は `monotonicity_condition_failed`、UE は U が `null`）、SP（1.0 eV は端の拒否で添字 [0]、2.5 eV は静的符号の `ValueError`）は表のとおりだった（135 ケース、想定外 0）。Codex も独立に、材料・エネルギー・指定・配列で束ねた 105 条件を確かめた（受理 88、端 6、定義域 9、静的符号 2）。
- **N2:** 単調性の条件 \((E^2-g^2)^3=125>(4/\pi)\Gamma gE^4\approx41.3\)、静的な \(F<0\) で、4 eV の事前検査は両指定で受理された。
- **固定深さの既存の挙動:** RJ4 の文言（LIM1 の `integration limit too small` と `loss_eV=5, index=(0,)`、\(\theta\) の `theta_deg`）は、固定深さの関数の現在の文言と同じ規則である。N1 の 10 Å・4 eV は固定深さで D3 の悪条件で止まり、N1 の 1 eV、P1 と P0 の 4 eV は成功した（Codex の確認）。
- **CLI の単点:** NC2（5 eV、IMP）と RJ3（1.999 eV）の設定は、現在の設定の読み込みで表の単点と精度を再現する（Codex の確認）。0.2–99.8 eV の格子を使うと、RJ3 の P1 は GC[9]（1.9999999999999998 eV）で端の拒否になり、別の理由で止まるので使わない。
- **Codex の確認で直したこと:** CLI の単点の設定（blocker 1）、V14(b) の入口ごとの入力と履歴の扱い（blocker 2）、V11(v) のための参照キー 44 個の追加と照合先の列挙（blocker 3）、EI の許す理由を族ごとにし、既存の失敗と理由の対応と優先順を定めた（major 4）、V7 を3値の判定にした（major 5）、V3(d) と V5 の単位と比の分数（major 6）、V1 の厳密な有理数と \(e_n\) の単独の検査（major 7）、直接和の診断の等式、例外の属性、V10(i) の口（major 8）、R1・R2 の手順と経路の値の凍結（major 9）、定例と全件の対応表と変異の表（major 10）、生成の検査（minor 11）、V2 と V15 の具体化（minor 12）。
- **Codex の再確認（2回目）で直したこと:** 前回の12件は解消3件、一部9件、未解消0件と判定された。新しい指摘は blocker 2件、major 5件。V14(b) の版番号を両方の版の比較から外して各版と照合し、失敗の対照で型が持つ全属性を比べ、UE の成功点を加えた（blocker B1）。V11(iii)・(v) を `check_plan` の全組（配列の MS10、MS3、MS4、MS15 を含む）で行う（blocker B2）。段階の順を M2-2 の順序（端 → 運動量の区間 → 単調性と静的符号 → 各点の証明）にし、複数の添字を返す拒否と1点を報告する未収束を分けた（M1）。精細化で減りうる係数の誤差を床から外し、失敗の側の観測（`panels`、評価の回数、`inner`）を加えた（M2）。`weight_out_of_range` を指数の範囲外に限り、零交差の直後は値を返すことと範囲外の対照を V1 に加えた（M3）。CLIDEFAULT と同じ精度の API の対照 CLIPREC を加え、CLI のエラーの欄を型ごとにした（M4）。照合の予定を検査・実行・添字の組で表に持たせ、予定と結果の件数を照合する（M5）。あわせて、参照の十進の保存の幅、経路1と経路3の具体的な規則と生成スクリプトの凍結、変異の表の検出する点と列、生成の独立の検査、V2 の尺度、V15 の再計算の検査を書いた。
- **Codex の最終確認（3回目）で直したこと:** 2回目の B1・B2 は解消、M1–M5 は解消4件・一部3件と判定され、明文の不整合が4か所示された。いずれも Codex の文案に沿って直した。K2 の検出先を、変異で値の変わる成分（W13 と W16e の虚部、W16h の実部）にした。経路1の積分域を \([0,r]\)（\(r=40\)、変える規則で 50）とし、幾何の端と位相の分割をその内部に限り、位相の分割を損失ごとの減衰の深さ \(t_c(T)\) で決めるようにした（「\(k_{\parallel-}>b\)」という誤った記述を直した。\(b\) より大きいのは全運動量の下端である）。尾の包絡をモード係数の絶対値から作り、作れない点は保留する。拒否の添字の規則を型ごとにした。V15 の再計算の対象を定例 1255、全件 1951 の MS の実行と定めた。あわせて、失敗の側の観測を `_panel_integrals` の記録との照合にし、`diagnostics["effective_settings"]`、経路3の `kpar_nodes` の変更、参照のキーの欄の一致の検査を加えた。この最後の修正の後の文書は、Codex の確認を受けていない。
- **表の再現:** 修正後の表の JSON は、Python 3.12.10（NumPy 2.5.5）と Python 3.10（NumPy 1.24.0）で同じ SHA-256 になった。

## 8. 記録

- 全件の記録（§3）の結果、所要時間、扱える精度の下限の観察を [検証記録](VALIDATION.md) に書く。
- 不合格の項目は、切り分けの結果とともに記録し、許容幅は変えない。照合未確定は、件数とキーを記録する。

## 付録　参照値と経路の値の凍結

（R1 と R2 が、実装の結果と比べる前に記入する: コミット、`references_m3a.json` と `paths_m3a.json` の SHA-256、生成の環境、所要時間、点数、不確かさの最大、照合未確定の点、未取得の点。）

### R1　多倍長の参照（参照担当 R1、2026-11-25。実装のコードと結果を見る前に記入）

- `benchmarks/numerics/m3a/references_m3a.json` の SHA-256: `1343541f939eb6d731f5957be48e6392698d9f27ecd7903675d57cde8c50e6f0`
- **コミット:** 凍結のコミット `f79d800` の作業ツリー（ブランチ `m3a-ref`）で作った。参照値を含むコミットは統括担当が記入する。
- **生成:** `benchmarks/numerics/m3a/generate_references_m3a.py`（mpmath と標準ライブラリだけ。`xpsloss` を読み込まない。表の JSON の SHA-256 を確かめてから読む）。比較の試験は `tests/test_m3a_reference.py`。
- **環境:** Windows 11（10.0.26300）、Python 3.12.10、mpmath 1.4.1（backend python）。12 プロセス（10 と 2 の2つの実行で同じキャッシュを共有）。同じ機械で他の担当の計算が動いていた。
- **所要時間:** 壁時計 3547 秒（約 59 分）、CPU の合計 39364 秒（経路 (a) 11661 秒、経路 (b) 22785 秒、残りは V3 のバルクの積分と極限の量）。最も長い群は Si・`k_total`・1.0 eV（7つの長さと極限）の 962 秒。JSON は所要時間を含むので、再生成すると値は同じでも SHA-256 は変わる。
- **点数:** 参照のキー 498（材料・エネルギー・指定・損失の 164 群。ギャップ下 130、ギャップより上 368）、極限のキー 192、V1 の点 35（うち深さ求積で確かめた点 8、入力の箱の8隅を持つ点 5）、V3(a)(b) の独立な 13a は 104 群（解析 69、バルクの積分 35）。
- **不確かさ:** \(u_j/B^{\rm ref}_j\) の最大は \(1.9\times10^{-10}\)（Al、`k_total`、20 Å、0.05 eV、`total`、TIGHT、\(Q_h=1\)。EI2 の点）。全 498 キーで、そのキーを使う実行のうち最も厳しい精度の組と電荷で \(u_j\le0.01B^{\rm ref}_j\) を満たす。(a) と (b) の相対差の最大は \(6.5\times10^{-19}\)（\(\lambda=10^{-6}\)、\(10^{-5}\) Å の `term_14`）、1 Å 以上の長さでは \(2.3\times10^{-21}\) 以下（(b) の尾 \(e^{-52}\) で決まる）。群の別の組立てと8項の和の相対差の最大は \(6.2\times10^{-45}\)。ギャップ下の単位電荷の `total` の生の8項の残差の最大は \(4.5\times10^{-39}\) Å⁻¹ eV⁻¹。
- **合成の式:** 平均した項は \(u_j=|a-b|+2|A-B|+r_A+\max(|\Delta_{\rm 次数}|,|\Delta_{\rm 位相}|,|\Delta_{\rm 尾}|)+\) 十進への丸め。13a は \(2|A-B|+r_A+\) 丸め。群は「別の組立てと (b) の項の和の差」\(+2|A-B|+r_A+\) 構成する項の (b) の変化の和 \(+\) 丸め。\(A\)、\(B\) は 50 桁の適応の実行と、80 桁で全パネルを二分した実行、\(r_A\) は A の受理した残差。内訳はキーと列ごとに JSON の `u_parts` にある。
- **照合未確定の点:** なし。**計算できなかった点:** なし。
- **経路 (b) の作り方についての注記:** 全キーの (b) は、印刷どおりの波括弧の基本の深さ因子（\(1\)、\(e^{-i\Omega_0a-k a}\)、\(e^{-2ka}\)、\((e^{-i\Omega_0a}-e^{-ka})(e^{i\Omega_0a}-e^{-ka})\)）を各深さの節点で数値評価して求積し、深さによらない因子を後で掛けた（費用のため）。印刷どおりの因子の全体を各節点で評価して求積した値との一致を、4条件（W4・`k_parallel`・√135 eV・20 Å、Al・`k_total`・15.6 eV・20 Å、Si・`k_total`・0.2 eV・5 Å、S0・`k_parallel`・5 eV・80 Å）で確かめた（相対差 \(7\times10^{-50}\) 以下。JSON の `self_checks`）。
- **統括担当の判断（2026-11-25、実装の結果を見る前）:** この作り方を受け入れる。§5 の (b) の趣旨は、(a) の平均済みの重みを使わずに、深さ方向の数値求積で (a) の閉形式を独立に確かめることである。基本因子への組替えは深さによらない係数だけの線形の書き換えで、各深さで印刷形と相対 \(1.5	imes10^{-60}\) で一致し、4条件の全体の評価とも一致した。全キーを全体の評価でやり直すことはしない。

### R2（経路1と経路3）

記入: 2026-11-25、R2。実装の結果とはまだ比べていない（実装はこの時点で存在しない）。

- **コミット:** 凍結のコミット `f79d800` の上の作業ツリー（ブランチ `m3a-paths`、未コミット。統括担当がコミットする）。
- `benchmarks/numerics/m3a/paths_m3a.json` の SHA-256: `42647fabcd07511f2c115ebaef4fa56b9bef22810b9856d8bf7aa3d02b489677`（3206393 バイト、ASCII、LF）
- **生成:** `benchmarks/numerics/m3a/generate_paths_m3a.py`（SHA-256 `abf8878ebd7695d80b07bb98eb94c92699d828c6e85f8ca5c02ba6f894e36b46`）。規則の実装は `tests/m3a_paths.py`（SHA-256 `e0c7fc468c043e1c09525f6c807296f12c2ce188ad19e1c29b5141a5c60da8a3`）。経路3は凍結した `tests/xps_route_b_reference.py`（変更なし）。
- **環境:** Python 3.12.10、NumPy 2.5.5、SciPy 1.18.1、xpsloss 0.2.0.dev0（`src` は `f79d800`）、Windows 11、Intel Core Ultra 7 155U（他の担当と共有、同時に最大10プロセス）。
- **所要時間:** 経路1の節点の評価の合計 182424 秒（Al 1000 eV 34113、Y0 13186、Si 133879、Al 4000 eV 1247）、予備の計算 1166 秒、尾の包絡 69 秒。経過時間はおよそ 5.7 時間（08:27–14:08）。
- **点数（経路1、V8 の全集合）:** 56 組、656 点。経路の基準（全節点・全列の \(\alpha_ie_{ij}\le f_i\eta D_j\) と \(u_j\le D_j/3\)）を満たした点 630、節点の失敗による経路の失敗 26、保留 0、未取得 0。
- **必須の3点:** すべて成功。\(\max_ju_j/D_j\) は (80 Å、99.8 eV) で \(4.90\times10^{-4}\)、(20 Å、1.0 eV) で \(9.61\times10^{-4}\)、(20 Å、15.6 eV) で \(4.00\times10^{-4}\)。
- **不確かさの最大:** 成功した点の \(\max u_j/D_j=6.12\times10^{-3}\)（Al、`k_total`、\(Q_h=1\)、80 Å、0.2 eV）。配分の比 \(\alpha_ie_{ij}/(f_i\eta D_j)\) の最大 0.896。\(\max(|\Delta_A|,|\Delta_B|,|\Delta_C|)/D_j\le1.3\times10^{-8}\)、尾 \(/D_j\le1.1\times10^{-6}\)。
- **経路の失敗（26 点、V8 は未完了）:** 0.05 eV の Al（\(Q_h=1\)、両指定、3長さの6点）と Y0（16点。\(Q_h=1\) の 80 Å の2組を除く全組）は、固定深さの「200 パネルで収束しない」で節点が失敗した。Si の 0.2 eV（\(Q_h=1\)、両指定、5 Å と 20 Å の4点）は「hard-gap arithmetic error exceeds an output budget」で節点が失敗した。いずれも節点に配分した精度が固定深さの下限を下回ったためで、点の集合からは外していない。尺度の最大を取る点が失敗した組では、その点だけ予備の計算の値で尺度を作った（JSON の `scale_source`）。
- **経路3（V9）:** 未取得（全 108 点、`not_computed`）。この計算機での計算量の見込みが 20 Å だけで約6日、全件で約9日（10並列）になり、実行していない。生成スクリプトは再開できる形で計画を持つ。完成させた場合は新しい SHA-256 をここに記録してから実装と比べる。
