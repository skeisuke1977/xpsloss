# M4b-2 の事前登録　束縛エネルギーの入出力と測定の記録

作成日: 2026-11-28。状態: **凍結（利用者の了解、2026-11-28。第4稿を凍結した）。**実装も参照値もない段階で固定した。Codex（実装担当）が第1稿（判定「修正後に凍結」。blocker 2、major 6、minor 1）、第2稿（同じ判定。新しく blocker 1、major 4、minor 2）、第3稿（同じ判定。blocker 0）を確かめた。この第4稿は、3回目の「残る指摘」だけを反映した（§10）。利用者が了解し（2026-11-28）、統括担当が凍結した。設計 §5 第6版の BE の台の判定（拒否せず弱い側の主張を返す）と、§9.1 の M4b-1 の試験2か所の変更も、利用者が凍結の前に了解した。この文書は、実装と参照値の計算の前に、検証の条件、許容幅、期待する結果と、実装が守る面（API、CLI、記録の形）を固定する契約になる。

**根拠:**

- [M4b の設計](M4B_DESIGN.md)（**第6版**）の §7（束縛エネルギーの入出力）、§5 の BE の台の判定（第6版で「食い違えば拒否」から「両方の側で求めて主張を弱い側に取る」に改めた）、§9 の `energy_axis`・例外・警告・CLI・`run.json`、§10 の「BE 軸」の群と変異 X25–X27・X32、§12 の U6・U9。
- [M4b-1 の事前登録](M4B1_PREREG.md)（`db67d12` で凍結）と、その表 `benchmarks/numerics/m4b1/expected_states_m4b1.json`（SHA-256 `d8d0f892…35e2`）。M4b-1 は §1 で BE を M4b-2 に回し、§4.1 で `energy_axis` を置かなかった。
- 手本は [M4a-2 の事前登録](M4A2_PREREG.md)（新しい数値の層がない段。基準のコードの子プロセスとのビット一致）と M4b-1 の事前登録（入力の書き方、「値なし」、変異の数え方、付録の自己参照）。
- Codex の確認 `work/consult/m4b2/prereg/codex_prereg_review.md`（1回目）、`codex_prereg_review2.md`（2回目）、`codex_prereg_review3.md`（3回目）。型と欄の名前、台の判定の規則、上界の比べ方、変異の経路、CI の行の選び方は、この確認の案で決めた。各指摘への対応と残る未決は、非公開の作業フォルダーの `work/consult/m4b2/prereg/open_questions.md` にある。

**範囲:** 束縛エネルギー（BE）の入出力（`multiple_loss_spectrum` の `energy_axis` と、CLI `xps-multiple-loss` の新しいスキーマ）と、測定の記録の形。内部と既存の API は運動エネルギー（KE）のまま（設計 §7）。`multiple_loss_order_areas`、M4a・M4a-2 の入口には足さない（`max_order=1` の1次は M4a とビット単位で同じなので、1次だけの BE もこの入口で得られる）。実測との比較、準位ごとの実行の合算、帯の比、基線は M4b-3。法線出射だけ。QUEELS-XPS の完全モデルは対象でない。

**重さ:** M4a-2 と同じ型。期待する状態の表、実装と別の担当が作る参照（R-C。パッケージを読み込まず、有理数で KE の値、ビンの所属、両側の台を固定する）、必須の変異。新しい数値の層がないので、試走と実材料の独立な参照（M4b-1 の R-B に当たるもの）は置かない。**物理の値の正しさは、BE の実行と、参照が独立に換算した KE の軸での同じ入口の実行とのビット一致で確かめる**（§7.1 の「物理量の対照」）。打ち切りの欄は、両側の規則（§2.3）で別に確かめる（「打ち切りの対照」）。KE の実行の正しさは M4b-1 の検証に依る。

**凍結の後の規則:**

- この文書と表を含むコミットの後に、参照値の計算と実装を始める。
- 比較の後に、許容幅と期待する結果を変えない。変える場合は、利用者の了解を得て、旧版を残し、別の試験として記録する（§11）。
- 必須の検査が時間の上限（§8）に収まらなければ、§8 の順で全件に回す。それでも収まらなければ旧版を残して再登録する。
- 不合格なら、(1) `energy_axis=None` の経路と既存の入口（基準のコードとのビット一致）、(2) BE の入口の検査（構造、値、排他）、(3) 換算と丸め、(4) 単射、(5) 両側の台と上界の幅、(6) KE での計算（M4b-1 の層。物理量の対照が通るか）、(7) 出力の向きと外の名前の復元、(8) 記録・CLI・`run.json`、(9) 参照の経路、の順に切り分けて記録する。

記号は設計の「表記」と M4b-1 の事前登録に従う。\(B\) は BE（フェルミ準位を基準に正の向き）、\(E\) は KE、\(C\) は換算の定数（\(E=C-B\)）、\(h\nu\) は光子エネルギー、\(\phi\) は `kinetic_energy_offset_eV`、\((B_\star,E_\star)\) は対応点、\(B_{\rm ref}\) は一次線の代表位置。\({\rm RN}(x)\) は実数 \(x\) を最も近い binary64 に1回丸めた値（偶数への丸め）、\({\rm Q}(x)\) は binary64 の \(x\) を厳密な有理数として読んだ値、\(\lceil x\rceil_{64}\) は \(x\) 以上の最小の binary64。\(\epsilon=2^{-52}\)。

## 1. 入力

数値は十進で与え、binary64 に変換して渡す。機械で読める同じ内容（十進と16進）が `benchmarks/numerics/m4b2/expected_states_m4b2.json` の `axes`、`records`、`binding_primaries`、`binding_output_edges`、`kernels` にある。**測定に見える値（線源、校正の線、準位）はすべて合成（synthetic）である。**文献の値や実測の値として使わない。

### 1.1 材料、解析の核、KE の入力

- **材料:** Al は `benchmarks/inputs/al_elfdb_v1.0.json`（SHA-256 `17ded0465a5e42124eebf44d9573369a5cbea87df59b5df11216ea2f7432d1c4`）。実行用の設定は `al_elfdb_v1.0.bulk.json`（`978931c1651881e6c7dfa739f4d626e4496225bf24bceb815466ff651934eae2`）。表の生成は両方のハッシュを照合する。`synthetic_wide`（\(g=0\)、\((20,200,8,0)\)、synthetic）は M4b-1 と同じ。
- **解析の核:** M4b-1 の K_A1s（窓 \([2,5]\)、13a \(=0.1\)）、K_A3（窓 \([1,8]\)、4つの箱。`hole` に負のビン）、K_T1（窓 \([1,2.5]\)、13a \(=0.2\)）、K_T2（窓 \([1,50]\)、13a \(=0.0072\)、15-T1 \(=0.005\)）を、凍結した M4b-1 の表からそのまま写す。新しい核は K_T1e（K_T1 の窓の下端を \(1+2^{-45}\) eV にしたもの）だけ。
- **KE の入力（基準のコードとの比較用）:** M4b-1 の D1000、E_M、E_A3 と、CI 用の小さい端 E_Ms（995、999、1000、1001 eV）。
- **共通の条件:** 公開の入口は \(L_w=20\) Å、`k_total`、窓 \([1,50]\) eV、乗数は `weight_length`、出射角 0°、\(E_0\) は行ごと。解析の入口は \(L_w=m=2\) Å、\(Q_h=1\)。書いていない引数は M4b-1 の既定。

### 1.2 軸（`energy_axis` の数値。すべて synthetic）

| 記号 | 形 | 与える値（eV） | \(C\)（eV） | 使い道 |
| --- | --- | --- | --- | --- |
| AX_PO | `photon_and_offset` | \(h\nu=1486.6\)、\(\phi=4.5\) | \({\rm Q}(1486.6)-9/2\)（binary64 で表せる） | 基本の行 |
| AX_PA | `paired_anchor` | \(h\nu=1486.6\)、\(B_\star=84.0\)、\(E_\star=1398.25\) | \(1482.25\)。含意されるずれは \({\rm Q}(1486.6)-1482.25\) | 対応点（AX_PO と別の \(C\)） |
| AX_R | `photon_and_offset` | \(h\nu=1486.6\)、\(\phi=4.6\) | \({\rm Q}(1486.6)-{\rm Q}(4.6)\)（**binary64 で表せない**。\({\rm RN}(C)=1482\)） | 2回の丸めの検出 |
| AX_D | `photon_and_offset` | \(h\nu=1253.5\)、\(\phi=4.25\) | \(1249.25\) | KE の像が M4b-1 の入力になる行 |
| AX_DA | `paired_anchor` | \(h\nu=1253.5\)、\(B_\star=249.25\)、\(E_\star=1000\) | \(1249.25\)。含意されるずれ 4.25 | CLI |
| AX_S | `photon_and_offset` | \(h\nu=130.5\)、\(\phi=2.5\) | \(128\) | 1 ulp の端が厳密に写る行 |
| AX_T2 | `photon_and_offset` | \(h\nu=1069\)、\(\phi=4.5\) | \(1064.5\) | KE の像が M4b-1 の T2c になる行 |
| AX_H | `photon_and_offset` | \(h\nu=6000\)、\(\phi=0\) | \(6000\) | 合成の HAXPES（設計 §8、U9） |
| AX_8192 | `photon_and_offset` | \(h\nu=8192\)、\(\phi=0\) | \(8192\) | 設計 §5・§7 の反例 |
| AX_16 | `photon_and_offset` | \(h\nu=16\)、\(\phi=0\) | \(16\) | 反例の端が厳密に写る対照 |

**記録（6つの欄。すべて synthetic）:** REC_FULL（全部を書いた記録。校正のずれ 0、帯電の補正なし、`emission_angle_deg` 0、`acceptance_half_angle_deg` 5）、REC_NR（6つの欄がすべて `"not_recorded"`）、REC_TILT（出射角 30）、REC_NAMED（線源の名前が `"Al Kalpha monochromatic"`）、REC_CORR（校正のずれ −0.25 eV、帯電の補正 1.5 eV を記録）、REC_NEST（`source.monochromatic` だけ `"not_recorded"`）。

### 1.3 BE の一次線と出力の端（eV。強度は \(I\)。\(u=2^{-40}\) eV）

| 記号 | 内容 |
| --- | --- |
| BD_2、BD_2r | δ 72.5（強度 1）と 60.5（強度 0.5）／並びを逆にしたもの |
| BB_neg、BB_neg_inc | ビン端 1.5、0.5、0、−0.5、−1.5（減少の順、負の BE、不等間隔）、積分強度 0.1、0.4、0.3、0.2／同じ一次線を増加の順で（強度も一緒に逆） |
| BD_ends、BD_dup | δ 70、74、76、強度 0.25、0.5、0.25／δ 74、74、70、強度 0.25、0.25、0.5（同じ位置の2本） |
| BD_731、BD_723、BD_249、BD_645、BD_0、BD_m1 | δ 73.1／72.3／249.25／64.5／0／−1（強度 1） |
| BD_249zero、BD_eps、BD_eps0 | δ 249.25（強度 0）／δ \(1+2^{-45}\)（強度 1）／同じ位置で強度 0 |
| BB_u | ビン端 0、\(2u\)、積分強度 1（設計 §5 第6版の「丸めた KE の問題」の例） |

| 記号 | 出力の BE の端 |
| --- | --- |
| BE_E1、BE_E1r | 70、72.5、73、74.5、75.5、77.5、80、81（不等間隔。7 ビン）／減少の順 |
| BE_Eneg | 14、10、7、4、2.5、1、0.25、−0.75（減少の順、負の端） |
| BE_E4 | 70、72、74、76 |
| BE_E5a、BE_E5、BE_E8 | 70.3、73.1、75.2、77.7、80.9、85.4／113.7、93.4、89.9、85.2、78.3、74.6、73.1、72.4（減少の順）／70.1、72.3、74.6、77.8、80.3、85.6（十進。換算で丸めが起きる） |
| BE_EM、BE_EMr、BE_EMs | 248.25、249.25、250.25、254.25、261.25、265.25、269.25、289.25（\(1249.25-\)E_M）／減少の順／248.25、249.25、250.25、254.25（\(1249.25-\)E_Ms。CI 用の小さい端） |
| BE_B120 | \(1249.25-(880.5+i)\)、\(i=0..120\)（減少の順。KE の像は M4b-1 の E_B120） |
| BE_T1_eq、BE_T1_dn、BE_T1_up | 63.5、64.5、65.5、66.5、67.5／最後の点を `nextafter(67.5, -inf)`／`nextafter(67.5, +inf)` |
| BE_T2 | 63.5、64.5、65.5、74.5、84.5、94.5、114.5（\(1064.5-\)E_T2） |
| BE_c7 | 0、1、2（設計 §7 の反例） |
| BE_c5、BE_c5b、BE_c5rev、BE_c5q | 0、\(2+2^{-43}\)／0、\(2+2^{-43}\)、9／0、\(2-2^{-43}\)／0、\(2.25+2^{-43}\) |
| BE_u | 0、\(3u/4\)、\(2u\) |
| BE_grid | `{"min": -2, "max": 2, "step": 0.5}`（CLI の短い形。負の `min`） |
| BE_neg_ke、BE_bad | 126、130（AX_S で KE が −2 eV）／70、72、71（単調でない） |

## 2. 凍結する規則

### 2.1 換算と丸め

- \(C\) は入力の binary64 から厳密な有理数で作る。`photon_and_offset` は \(C={\rm Q}(h\nu)-{\rm Q}(\phi)\)、`paired_anchor` は \(C={\rm Q}(B_\star)+{\rm Q}(E_\star)\)。**\(C\) を binary64 に丸めてから使わない。**
- 各座標（一次線の位置・ビンの端、出力の端）の KE は \({\rm RN}(C-{\rm Q}(B))\)。入力の数から \(C-B\) までは有理数で、丸めは最後の1回だけ。
- **返す数値は、丸めた KE の座標の問題の値である**（設計 §5 第6版）。換算の丸めは返す誤差に入れない。δ の所属は単射（§2.2）で保たれる。**ビンの一次線の積分強度は、換算した KE の幅に一様に配る。**したがって、丸めが起きるビンの一次線では、ゼロ損失や各次数のビンの値が、元の BE の厳密な重なりから作った値と丸めの水準を超えて違いうる（§7.2 の例）。
- 記録の値: `conversion_constant_eV` は \({\rm RN}(C)\)（計算には使わない）、`conversion_constant_exact` は \(C\) の既約分数の文字列（「分子/正の分母」。整数も `p/1`）。`rounding_max_eV` は全座標の \(|{\rm RN}(C-B)-(C-B)|\) の最大を上へ丸めた binary64（すべて厳密なら `0.0`）。有理数が 0 の記録は `+0.0`。
- `paired_anchor` の含意されるずれ `implied_offset_eV` は \({\rm RN}({\rm Q}(h\nu)-C)\)。**記録するだけで、\(C\) に使わない。**`photon_and_offset` では `None`。
- `primary_reference_kinetic_eV` は \({\rm RN}(C-{\rm Q}(B_{\rm ref}))\)。`electron_energy_minus_reference_eV` は \({\rm RN}({\rm Q}(E_0)-(C-{\rm Q}(B_{\rm ref})))\)（解析の入口では `None`）。\(E_0\) は BE から決めない（`electron_energy_eV` をそのまま核に渡す）。
- 座標の KE が負か、有限でなければ `invalid_energy_axis`（`negative_kinetic_energy`、`not_finite`）。記録の値（\({\rm RN}(C)\)、含意されるずれ、代表位置の KE、\(E_0\) との差）が binary64 で表せなければ、その欄の `not_finite`。代表位置の KE と対応点は記録なので、負でも受け付ける。

### 2.2 単射の検査

- 対象は、**一次線の位置（強度 0 の線も含む）・一次線のビンの端・出力の端の全部を合わせた、異なる BE の座標の集合** \(X\)（binary64 として等しい値は1つに合併する。元から等しい座標と、同じ位置の δ の重複は許す）。\(B_{\rm ref}\) と対応点は入れない。
- \(X\) を昇順に並べ、隣り合う \(x<y\) で \({\rm RN}(C-x)>{\rm RN}(C-y)\) を求める。成り立たない組があれば `axis_conversion_not_injective`。**BE の昇順で最初の隣り合う組**を報告する。各座標の代表の出所（`roles`）は「一次線の入力の順、次に出力の端の入力の順」で最初に現れたもの（例 `"primary.binding_energies_eV[0]"`、`"output_edges_eV[1]"`）。
- 配列ごとの検査では足りない（設計 §7 の反例）。

### 2.3 両側の台と上界の幅（設計 §5 第6版。拒否しない）

Codex の確認の「採る規則」による。

1. **台の点と幅:** \(P\) は正の強度の δ の位置と、強度が正の一次線のビンの両端の、異なる点の集合。\(P\ne\varnothing\) のとき、表示のビン \(i\)（入力の順）の高い BE の端を \(d_i\) として、\(U_{B,i}={\rm Q}(B_{d_i})-\min_P{\rm Q}(B_p)\)、\(U_{K,i}=\max_P{\rm Q}(E_p)-{\rm Q}(E_{d_i})\)（\(E\) は丸めた後の binary64）、\(U_i=\max(U_{B,i},U_{K,i})\)。
2. **台の次数:** 各側で \(S_j=\max(0,\lfloor\max_iU_{j,i}/{\rm Q}(a)\rfloor)\)（\(j\) は B、K）、\(S=\max(S_B,S_K)\)。\(P\) が空なら両側とも 0。`truncation.support_order` は \(S\)。整数の `max_order` は \(N\) を変えない。`"display_support"` は \(N=S\)（0 もある）。\(S\) は M4a・模型・設定の検査の後に決め、続いて `order_limit` と作業量の検査を行う（次数についてのループを上限の検査の前に置かない）。
3. **`complete` の条件:** 次数 \(k\) まで計算したとき、ビン \(i\) の上界が台を根拠に `(0.0, "complete", None)` になるのは、\(P\) が空か、**両方の側で** \(U_{j,i}<(k+1){\rm Q}(a)\) のときだけ。**等号は届く側。**比較と床は有理数のまま行う。
4. **上界の幅:** 届くビンの上界（M4b-1 §3.3 の式）の \(U_D\) には \(\lceil U_i\rceil_{64}\) を使う。**全部のビンで常に大きい方を取る**（食い違いの件数が 0 でも省かない）。上界の候補（\(\tau\) の格子、許されるときの \(\tau\to0\)）、外向きの丸め、最小と同点の規則、厳密なゼロの因子による `(0.0, None, τ)`、下へのあふれ、`no_usable_tau`、`not_representable` は M4b-1 のまま。上界は、そのビンから低い BE 側の全部に落ちる量の累積の上界である。
5. **状態:** \(N\ge S\) なら `complete_on_display`。そうでなく許容がなければ `truncated`。許容があれば、表示の全ビンの上界が有限で、最大が許容 \(\times S_{\rm abs}\) 以下（等号可）なら `bounded_on_display`、欠けた上界があれば `truncation_bound_unavailable`、ほかの超過は `truncation_bound_exceeds_tolerance`。`required_order_estimate` も \([0,S]\) で同じ両側の規則を使う。`neglected_total_area` は軸に依らない全軸の量で、KE の実行の値のまま。
6. **診断:** 各 \(p\in P\)、出力の端 \(d\)、\(1\le n\le N\)、\(c\in\{n{\rm Q}(a),n{\rm Q}(b)\}\) について、\({\rm Q}(B_d)-{\rm Q}(B_p)\) と \({\rm Q}(E_p)-{\rm Q}(E_d)\) の \(c\) に対する大小・等号（3通り）を比べ、違う組を数える（計算した次数の食い違い）。次の次数については、各ビンの \(U_{B,i}\) と \(U_{K,i}\) の \((N+1){\rm Q}(a)\) に対する大小・等号を比べ、違うビンを数え（次の次数の食い違い。等号と「大」の違いのように、両側とも届く場合も数える）、届く・届かないが違うビンを別に記録する。\((N+1)b\) は使わない。
7. **記録と警告:** `energy_axis` に `support_order_binding`（\(S_B\)）、`support_order_kinetic`（\(S_K\)）、`region_loss_upper_eV`（\(\lceil U_i\rceil_{64}\) の、入力の順の配列。\(P\) が空なら全部 `None`）、`support_comparison`（`computed_relation_mismatch_count`、`next_relation_mismatch_count`、`next_reachability_mismatch_bins`（入力の 0 始まりの添字）、`first_mismatches`）を置く。件数は全数。**`first_mismatches` は、計算した次数の食い違いと、次の次数の食い違いの両方を入れた先頭 8 件**で、順は「次数 → 境（`lower`、`upper`）→ 一次線の BE の昇順 → 端の BE の昇順」（次の次数の例は次数 \(N+1\)、境 `lower` なので、計算した次数の例の後に来る）。各例は、比較の種類（`computed_order`、`next_order`）、次数、境、一次線の点と端のそれぞれの BE と KE、元の入力の添字（`primary_index`、`edge_index`）、2つの関係を持つ。次の次数の例は、一次線の点に台の最低の BE の点、端にそのビンの高い BE の端を使い、入力のビンの添字 `bin_index` も持つ。\(S_B\ne S_K\) か、どれかの件数が 0 でなければ、警告 `axis_conversion_support_disagreement` を1件だけ出す。
8. **計算量:** 単射は並べ替えを含めて \(O(M\log M)\)（\(M=|X|\)）、台の次数と各ビンの幅は \(O(|P|+|D|)\)。診断は、各 \(p,n,c\) で BE の昇順の端に対する関係が両側とも単調なので二分探索で数えられ、\(O(N|P|\log|D|)\)、メモリ \(O(M)\)。製品は二分探索で件数を出し、例を全件は保存しない（表の生成と参照は、独立な全探索で数えてよい）。

### 2.4 出力の向き、外の名前、境の上の δ

- 計算は KE の昇順で行う。一次線のビンは KE の昇順に並べ直し、積分強度も一緒に逆にする。δ は重複を保ち、KE の昇順に安定に並べ、強度も同じ置換にする。入力の順は別に保存する。
- **出力は入力の BE の端の向きに戻す。**`bin[i]` は常に入力の端の \(i\) 番目と \(i+1\) 番目の間のビン。入力が BE の増加の順なら内部の KE のビン \(n_b-1-i\)、減少の順なら内部の KE のビン \(i\)。
- 外の領域は `outside.high_binding`（出力の最高の BE より高い側。内部の `outside.low`）と `outside.low_binding`（最低の BE より低い側。内部の `outside.high`）。入力の向きで名前は変わらない。
- **境の上の δ（出力の端を昇順に \(\beta_0<\dots<\beta_{n_b}\) として）:** 最低の BE のビンは \([\beta_0,\beta_1]\)、ほかのビンは \((\beta_j,\beta_{j+1}]\)、`outside.low_binding` は \((-\infty,\beta_0)\)、`outside.high_binding` は \((\beta_{n_b},\infty)\)。内部の境の上の δ は低い BE 側のビンへ、全体の両端の上の δ は端のビンへ入る。現行の `_fold.Geometry.mass` の規則（KE のビンは \([e_i,e_{i+1})\)、最高の KE のビンだけ右も閉）を \(E=C-B\) で写したもので、設計 §7 の文と一致する。
- **内部の KE の索引のまま残るもの:** `single_loss`（§4.3）、解析の入口の2つ目の戻り値 `panels_by_order`（損失の順）、`diagnostics` の中の領域の名前（`balance[n].components`、網の履歴など）、M4a・M4b-1 から通った例外の `key`・`exceeded` は、KE の `bin[i]`、`outside.low`、`outside.high` で書かれたままにする。BE に直すのは、§4.3 に挙げた欄と、BE の実行の打ち切りの例外が持つビンごとのリストだけ。

## 3. 測定の記録の形

**スキーマ `xpsloss-measured-spectrum-v0`**（CLI `xps-multiple-loss` の `--spectrum` のファイル。JSON、UTF-8）。API の `BindingEnergyAxis` と同じ欄を持つ。

| 欄 | 必須 | 型と単位 | `"not_recorded"` |
| --- | --- | --- | --- |
| `schema`、`axis` | 必須 | `"xpsloss-measured-spectrum-v0"`、`"binding_energy_eV"` だけ | — |
| `measurement.photon_energy_eV` | 必須 | 有限の正の数（eV） | **不可** |
| `measurement.reference.kind` | 必須 | `"photon_and_offset"` か `"paired_anchor"` | 不可 |
| `measurement.reference.kinetic_energy_offset_eV` | `photon_and_offset` で必須 | 有限の数（eV） | **不可** |
| `measurement.reference.binding_anchor_eV`、`kinetic_anchor_eV` | `paired_anchor` で必須 | 有限の数（eV） | **不可** |
| `measurement.reference.line`、`citation` | `paired_anchor` で必須 | 文字列（線の名前、出典） | 可（警告） |
| `measurement.primary_reference_binding_eV` | 必須 | 有限の数（eV） | 不可 |
| `measurement.primary_reference_choice` | 必須 | 文字列（選び方） | 可（警告） |
| `measurement.source` | 必須 | `{"name": 文字列, "monochromatic": 真偽, "satellite_treatment": 文字列}` | 欄の全体か、中の各値（警告） |
| `measurement.calibration` | 必須 | `{"line": 文字列, "reference_binding_eV": 数, "measured_binding_eV": 数, "applied_shift_eV": 数}`（eV） | 同上 |
| `measurement.charge_correction` | 必須 | `{"applied": 真偽, "method": 文字列, "shift_eV": 数}` | 同上 |
| `measurement.core_level` | 必須 | 文字列 | 可（警告） |
| `measurement.emission_angle_deg` | 必須 | 有限の数（表面の法線からの角度、度） | 可（警告） |
| `measurement.acceptance_half_angle_deg` | 必須 | 有限の非負の数（度） | 可（警告） |
| `primary` | 必須 | `{"kind": "delta", "binding_energies_eV": [...], "intensities": [...]}` か `{"kind": "binned", "binding_edges_eV": [...], "intensities": [...]}`。任意の欄は KE のスキーマと同じ5つ | — |
| `output_binding_edges_eV` | 必須 | 配列（狭義の増加か減少。負を許す）か `{"min", "max", "step"}` | — |

- **値の規則:** 文字列は UTF-8 で、空白だけのものは不可（`empty`）。数に真偽値を混ぜない（`type`）。\(h\nu\le0\) と負の取り込み半角は `out_of_range`。`charge_correction.applied` が偽なら `method` と `shift_eV` は `null`、真なら所定の型か `"not_recorded"`、`applied` が `"not_recorded"` なら2つとも `null`。記録の3つの辞書は、所定のキーだけを持つ。受け付けた記録は深く写して持つ（渡した辞書を後で変えても結果は変わらない）。
- **排他:** `reference` は、選んだ形の欄だけを持つ。どちらの形でも `photon_energy_eV` は必須で、線源の名前から数値を補わない（D10）。
- **`"not_recorded"` の一覧:** 結果の `energy_axis["not_recorded"]` は、欄の名前のタプルで、順は `reference.line`、`reference.citation`、`primary_reference_choice`、続いて記録の6欄を上の表の順、各欄の中は表の中の順（`source.monochromatic` の形）。親が `"not_recorded"` なら親の1件だけ。
- **記録は計算に使わない。**`calibration.applied_shift_eV`、`charge_correction.shift_eV`、`emission_angle_deg` は「何をしたか」の記録で、パッケージは軸をずらさず、角度も計算に使わない（設計 §7、U6。BW3）。真空準位を基準にした BE は扱わない。
- **観測の計数は入れない。**このスキーマは、一次線（模型への入力）、出力の端、測定の記録だけを持つ。実測のスペクトルそのもの（計数、誤差、透過、時間）の欄は置かず、未知のキーとして拒否する（BC05）。比較の入力は M4b-3 の事前登録で、データを見る前に決める。設計 §7 の「誤差の配列」は、この段では置かない。
- **短い形の端:** `{"min", "max", "step"}` は、有限、`min < max`、`step > 0`、`min` は負でもよい。既存の `_edges` と同じ手順（十進の文字列で区間の数が整数であることとセルの数の上限を割り付けの前に確かめ、`numpy.linspace(min, max, 区間の数 + 1)` で展開する）。増加の順になる。形の誤りは `invalid_schema`、展開した端の値の誤りは `invalid_output_edges`、資源の拒否は既存の `cell_limit_exceeded`。
- **判定の順（契約）:** (1) JSON の構造、(2) 軸と記録の値の意味、(3) 一次線と出力の端、(4) 換算と単射。CLI は (2) の前に BE の一次線の型を作らない。
  - (1) は `invalid_schema`: JSON として不正、重複したキー、最上位・`measurement`・`primary`・記録の辞書の未知のキーと必須のキーの欠落、未知の `schema`・`axis`・`primary.kind`、短い形の端の形の誤り。
  - `reference` のキーは2つの形の和集合を受け取る。**選ばなかった形のキーを書いていれば（値が `null` でも）`reference`・`exclusive`**、選んだ形の必須の値が `null` なら、その欄の `missing`。CLI は、キーの有無を dataclass の既定値に落とす前に調べる。API では、省略と既定の `None` は同じで、選ばなかった形の欄が `None` でなければ `exclusive`。
- **入力のハッシュ:** `run.json` の `spectrum_input` は `file`、`sha256`（ファイル全体）、`schema`、`binding_output_edges_sha256_float64_le`（入力の順の BE の端）、`kinetic_output_edges_sha256_float64_le`（内部の昇順の KE の端）を持つ。どちらも連続したリトルエンディアンの `float64` のバイト列の SHA-256。
- KE のスキーマ `xpsloss-single-loss-spectrum-v0` は変えない（`axis` は `"kinetic_energy_eV"` だけ）。理論の計算（KE を直接に与える）では、光源も BE の基準も求めない。

## 4. 実装が守る約束（試験が使う面）

実装担当は参照値を見ない。この節と §5 の表だけで実装と試験を書けるようにする。名前と形は Codex の確認の「凍結する名前と形の案」による。

### 4.1 公開の型と署名（`xpsloss.xps_multiple_loss` と `xpsloss` の両方から、同じオブジェクトを公開する）

```python
@dataclass(frozen=True, kw_only=True)
class EnergyAxisReference:
    kind: str                                   # "photon_and_offset" | "paired_anchor"
    kinetic_energy_offset_eV: float | None = None
    binding_anchor_eV: float | None = None
    kinetic_anchor_eV: float | None = None
    line: str | None = None
    citation: str | None = None

@dataclass(frozen=True, kw_only=True)
class BindingEnergyAxis:                        # 既定値なし（すべて必須）
    photon_energy_eV: float
    reference: EnergyAxisReference
    primary_reference_binding_eV: float
    primary_reference_choice: str
    source: Mapping | str
    calibration: Mapping | str
    charge_correction: Mapping | str
    core_level: str
    emission_angle_deg: float | str
    acceptance_half_angle_deg: float | str

BindingDeltaPrimary(binding_energies_eV, intensities, intensity_unit="arb. unit", description=None,
                    included_broadening=None, truncation=None, intrinsic_satellites=None)
BindingBinnedPrimary(binding_edges_eV, intensities, ...同じ任意の欄)

multiple_loss_spectrum(material, electron_energy_eV, primary, output_edges_eV, *, ...M4b-1 の引数...,
                       energy_axis=None) -> MultipleLossSpectrumResult
```

- **`energy_axis=None`（既定）は M4b-1 の経路そのもの**で、検査の順、演算の順、結果の全欄が基準のコードとビット単位で同じ（§5 の BASE の群）。旧い引数の位置、既定、戻り値の型を保ち、`energy_axis` は末尾のキーワード専用の引数。結果に足す欄は `energy_axis`（`None`）だけ。
- `energy_axis` が `BindingEnergyAxis` のとき、`primary` は BE の一次線の型、`output_edges_eV` は BE の端（狭義の増加か減少、2点以上、有限。負を許す）。
- **2つの軸の型は構築のときに検査しない**（値の検査は API の中で行い、理由を `invalid_energy_axis` にそろえる）。**BE の一次線の型は構築のときに検査する。**KE の型のサブクラスにはしない。位置・端・強度は、独立した1次元の `float64` の写しで、有限、読み取り専用。位置と端は符号を許し、端は狭義の増加か減少。強度の長さ・上下限・任意の文字列は現行の `_single_loss_types._primary` と同じ。違反は `SingleLossInputError("invalid_primary")`。
- `multiple_loss_order_areas` と M4a・M4a-2 の入口は `energy_axis` を取らない（渡せば `TypeError`）。既存の関数、型、例外、CSV、`run.json`、既定の数値は変えない。

### 4.2 検査の順と例外

API の順は、構築の済んだ引数についてのものである（BE の一次線の型の誤りは、構築のときに先に出る）。

| 順 | 検査 | 理由（型） |
| --- | --- | --- |
| 0 | `energy_axis is None` なら M4b-1 の経路（BE の一次線の型を渡せば、現行のとおり `invalid_primary`） | `SingleLossInputError` |
| 1 | `energy_axis` の型、数値、`reference` の形と排他、`"not_recorded"` の可否、記録の欄の形（§3） | `invalid_energy_axis`（`MultipleLossInputError`） |
| 2 | `primary` が BE の型であること | `invalid_energy_axis`（`field="primary"`） |
| 3 | 出力の BE の端の形（単調、有限、2点以上） | `invalid_output_edges`（`SingleLossInputError`。M4a の理由のまま） |
| 4 | 換算（§2.1）。負か非有限の KE、表せない記録の値 | `invalid_energy_axis` |
| 5 | 単射（§2.2） | `axis_conversion_not_injective` |
| 6 | M4a の入力の検査と M4b-1 の `invalid_model`、`invalid_max_order`、`invalid_neglect_tolerance`、`invalid_settings`（既存の `_run`／`on_validated` の順のまま） | M4a・M4b-1 のまま |
| 7 | \(S\) を両側で決める（§2.3）。続いて `order_limit_exceeded`、`work_lower_bound_exceeded`（M4b-1 のまま。次数は \(N\)、`"display_support"` では \(S\)） | M4b-1 のまま |

- 1–5 は KE の型を作る前（M4a の実行の前）に行う。解析の入口でも、`energy_axis` を与えたら、先頭の窓の検査（`checked_window`）より前に 1–5 を行う。**2つの新しい理由は、どちらも核の評価 0 回。**
- **台の食い違いは拒否しない**（第5版の `axis_conversion_changes_support` は公開の理由にしない）。BN19–BN22 は、食い違いのある入力でも M4b-1 の理由がそのまま出ることを見る。
- 1 の中で複数の欠陥が重なるときの `field` の優先は登録しない（各行の欠陥は1つ）。

**属性**（該当しないものは持たないか `None`。**新しい属性は BE の理由のときだけ付け、既存の理由の例外の `vars(exc)` と CLI の失敗の JSON は基準と同じにする**。BASE-n13）:

| 理由 | 属性 |
| --- | --- |
| `invalid_energy_axis` | `field`（点でつないだ欄の名前。`"energy_axis"`、`"photon_energy_eV"`、`"reference"`、`"reference.kind"`、`"source"`、`"primary"`、`"output_edges_eV"` など）、`detail`（`"type"`、`"missing"`、`"not_finite"`、`"out_of_range"`、`"empty"`、`"not_recorded_not_allowed"`、`"exclusive"`、`"unknown_kind"`、`"requires_binding_primary"`、`"negative_kinetic_energy"`） |
| `axis_conversion_not_injective` | `binding_eV`（衝突した2つの BE、昇順）、`kinetic_eV`（その2つの丸めた KE、同じ順）、`roles`（2つの座標の出所）。API では2要素のタプル、JSON では配列 |
| `invalid_schema`（CLI） | `detail`（M4b-1 は文言にだけ使っていた。新スキーマでは属性として公開する） |

### 4.3 結果

`MultipleLossSpectrumResult` に欄 `energy_axis` を足す（KE では `None`）。BE では:

- **入力の BE の向きの欄:** 最上位と `orders` の配列（`zero_loss`、`spectrum`、`<c>`、誤差、誤差の内訳）、`truncation` のビンごとの3つのリスト。
- **外の領域:** `outside` と `orders[n].outside` のキーは `"low_binding"`、`"high_binding"`（`"low"`、`"high"` は持たない）。**最上位の `outside[side]` は `zero_loss`、`zero_loss_error`、`spectrum`、`spectrum_error` の4つを持つ**（KE の結果の最上位の `outside[side]` は `spectrum` と誤差だけで、外のゼロ損失は `single_loss.outside[side]` にある。BE では内包した KE の結果を読ませない）。
- **表のキーの読み場所:** `zero_loss.bin[i]` は `result.zero_loss[i]`、`zero_loss.outside.<side>` は `result.outside[side]["zero_loss"]`、`order[n].bin[i].<c>` は `getattr(result.orders[n-1], c)[i]`、`order[n].outside.<side>.<c>` は `result.orders[n-1].outside[side][c]`。KE の対照の外のゼロ損失だけは `control.single_loss.outside[side]["zero_loss"]` を読む。
- **軸に依らない欄:** 面積、感度、\(\hat r\)、`mesh`、`neglected_total_area` は KE の実行と同じ。
- **`single_loss` は内部の KE の結果のまま**（KE の昇順、`outside` は `"low"`・`"high"`、`primary` は換算した KE の型、密度の分母は KE の幅）。「M4a の結果をそのまま」（設計 §9）を保つ。docstring に明記する。
- 解析の入口の戻り値 `(result, panels_by_order)` は、1つ目だけを BE に直す。
- `energy_axis` は辞書で、次のキーを持つ。

| キー | 内容 |
| --- | --- |
| `axis`、`binding_reference` | `"binding_energy_eV"`、`"fermi_level"` |
| `reference_kind`、`photon_energy_eV`、`kinetic_energy_offset_eV`、`binding_anchor_eV`、`kinetic_anchor_eV`、`anchor_line`、`anchor_citation` | 入力の写し（該当しない形の欄は `None`） |
| `conversion_constant_eV`、`conversion_constant_exact`、`implied_offset_eV`、`rounding_max_eV` | §2.1 |
| `primary_reference_binding_eV`、`primary_reference_choice`、`primary_reference_kinetic_eV`、`electron_energy_minus_reference_eV` | §2.1 |
| `input_binding_edges_eV`、`input_primary`、`direction` | 入力の BE の端（`float64`、入力の順とビットを保つ）、入力の BE の一次線の型そのもの、`{"output_edges": "increasing"\|"decreasing", "primary": "increasing"\|"decreasing"\|None}`（δ は `None`） |
| `kinetic_output_edges_eV`、`kinetic_primary` | 換算した KE の端（昇順）と KE の一次線（`single_loss` のものと同じ値）。表のキー `energy_axis.kinetic_primary.coordinates[k]` は、KE の昇順の位置かビンの端 |
| `support_order_binding`、`support_order_kinetic`、`region_loss_upper_eV`、`support_comparison` | §2.3 |
| `record`、`not_recorded` | 6つの記録の欄の深い写しと、`"not_recorded"` の欄の名前のタプル（§3） |

### 4.4 警告（出所 `xps_multiple_loss`。M4b-1 の識別子の後に、この順で足す）

| 識別子 | 条件 |
| --- | --- |
| `binding_axis_assumed_calibrated` | `energy_axis` を与えたとき常に（入力は校正・補正の済んだフェルミ準位基準の目盛りとみなし、パッケージは軸をずらさない） |
| `energy_axis_field_not_recorded` | `not_recorded` が空でない（1件。文言に欄の名前を並べる） |
| `recorded_emission_angle_differs_from_theta` | 記録した `emission_angle_deg` が数値で、計算の `theta_deg` と違う |
| `axis_conversion_support_disagreement` | §2.3 の 7（1件だけ） |

外側の `warnings`（文言）と `diagnostics["warnings"]`（`source`、`code`、`message`）の順は対応させる。M4a の `E0_outside_primary_support`（換算した KE の台で判定）と M4b-1 の `truncated`、`E0_above_5keV` などは、そのまま出る（`truncated` は §2.3 の状態による）。\(E_0\) と \(C-B_{\rm ref}\) の差には警告の境を置かない（診断に残すだけ。設計 §7）。

### 4.5 試験用の口

`xpsloss._loss_orders.integrate_orders(..., energy_axis=None)` に同じ末尾のキーワード専用の引数を足す。意味は公開の入口と同じ（\(E_0\) がないので `electron_energy_minus_reference_eV` は `None`、記録した出射角は 0 と比べる）。BE の行の多くはこの入口で作る（CI の時間を抑える。核の呼び出しの回数も数えられる）。ほかの試験用の口は足さない。変異はレビュー担当が本番のコードに入れる。

**上界の上限の照合（§7.2）で、照合器が読む場所（既存の非公開の関数だけ。公開の API は足さない）:** 照合器は、KE の対照の実行の間、`xpsloss._loss_orders` が呼ぶ `_order_bounds.truncation(single, panels, geometry, requested, n, support, ratio, tol, negative, diagnostics, tail_model=None)` を読み取り専用で包み、渡った実引数（1次のパネル、乗数の区間を持つ `single`、次数）を写す。解析の入口でも公開の入口でも同じ場所を通る。上限は、写した実引数のまま、`geometry` を「幅が R-C の \(\lceil U_i\rceil_{64}\) になる代わりの入れ物」（`points=(w,)`、`edges=(0.0, 0.0)`）に替え、`tol=None`、`tail_model` に「M4b-1 の既定の式を同じ順に評価し、KE の対照の最適の \(\tau_K\) のほかの候補を \(+\infty\) にする」関数を渡して、同じ関数をもう一度呼んで作る。**この関数の署名と `tail_model` の意味（`(taus, count, logk, logr)` を受けて候補の対数の裾を返す）を、M4b-2 の実装で変えない。**製品が BE の実行でこの関数に渡す幅は §2.3 の 4 のとおり。

### 4.6 CLI、CSV、`run.json`、`capabilities`

- `xpsloss xps-multiple-loss CONFIG --spectrum SPECTRUM.json ...`（引数は M4b-1 のまま）。`SPECTRUM.json` の `schema` で分ける。旧スキーマ `xpsloss-single-loss-spectrum-v0` は KE だけのまま（`axis` が `"binding_energy_eV"` なら、これまでどおり `invalid_schema`。**M4b-1 の C09 は変えずに済む**）。新スキーマは `xps-multiple-loss` だけが受け付ける（`xps-single-loss` は `SingleLossInputError("invalid_schema")`）。
- **理由の分け方と順**は §3 の「判定の順」。構造の誤りは `MultipleLossInputError("invalid_schema")`、値の意味は API の検査に渡して API の理由で返す。
- **失敗:** 終了 2、出力先を作らない、標準エラーに JSON 1つ（`status`、`error_type`、`reason`、`message` と §4.2 の属性。浮動小数は読み直して同じ binary64 になる桁数）。CLI の文字列は ASCII。
- **`multiple_loss_spectrum.csv`（BE の実行）:** 見出しは `region,binding_edge_low_eV,binding_edge_high_eV,kinetic_edge_low_eV,kinetic_edge_high_eV,zero_loss,zero_loss_error`、続いて M4b-1 と同じ次数ごとの列、最後に `spectrum,spectrum_error`。4つの端の列は、その領域の下限と上限を数値の順で書く。`low_binding` の行は BE の下限と KE の上限が空欄、`high_binding` の行は BE の上限と KE の下限が空欄。行の並びだけが入力の向きに従う（増加なら `low_binding`、`0`…、`high_binding`、減少なら `high_binding`、`0`…、`low_binding`）。**値は BE の結果の最上位と各次数から作る**（内包した M4a の配列をそのまま列にしない）。`%.17g`。`order_areas.csv` は M4b-1 のまま。出力のファイルの名前は3つのまま。
- **`run.json`（BE の実行）:** M4b-1 の欄に加えて、最上位の `energy_axis`（§4.3 の辞書を JSON にしたもの。`results.energy_axis` と同じ内容）、`measurement`（ファイルの `measurement` の写し）、`spectrum_input`（§3 の5つ）、`spectrum_definition.axis`（`"binding_energy_eV"`）、`units` の `binding_energy`（`"eV"`）。最上位の `primary` と `output_edges_eV` は入力の BE の記録で、`axis` を明示する。`results` は §4.3 の向きと名前。`results.single_loss` の `primary` と端は KE。`energy_axis.input_primary` は、`kind`、位置か端、強度、任意の5欄だけを持つ、そのまま再入力できる辞書（計算した積分や台を入れない）。`effective_settings` は 13 個のまま。`comparison` は常に `not_matched`、`physical_validity` は `not_established`。
- **KE の実行の `run.json`、CSV、失敗の JSON は変えない。**KE の `run.json` に `energy_axis` のキーをどこにも出さない（`results.energy_axis: null` も出さない。BASE-cli）。
- **再生（BC01s、BC01、BC02、BC07、BC12）:** 書いたファイルだけから入力を作り直す。使う欄は、`config`、最上位の `model_id`、`measurement`、`energy_axis.input_primary`、`energy_axis.input_binding_edges_eV`、`effective_settings`、`spectrum_definition`（窓、運動量、角度、正孔電荷、重みと乗数の選択、\(E_0\)、`max_order`。導出の長さでは API の `weight_length_A=None` と種・追加の逆数の長さ、外部の乗数では `multiplier_A` に対応づける）。`energy_axis` と `measurement` の重なる欄は、意味と値を照合してから使う。API を呼び、値、誤差、状態、警告を照合し、2つの CSV の全部のセルを独立な列の対応で照合する。BE と KE の2つの端のハッシュも独立に計算して照合する。
- **`capabilities`:** `multiple_scattering.unsupported_physics` から `binding_energy_axis` を外す（ほかの9個は順もそのまま）。同じブロックに `energy_axes`（`["kinetic_energy_eV", "binding_energy_eV"]`）と `measured_spectrum_schema`（`"xpsloss-measured-spectrum-v0"`）を足す。`supported` と、ほかのブロックは変えない。

## 5. 期待する状態（82 行）

**表:** `benchmarks/numerics/m4b2/expected_states_m4b2.json`（スキーマ `xpsloss-m4b2-expected-states-v1-draft4`、ASCII、LF。SHA-256 は凍結の時に統括担当が付録に書く）。生成は `python benchmarks/numerics/m4b2/expected_states_m4b2.py`（標準ライブラリだけ。`xpsloss` を読み込まない）。`--check` で再生成してバイト単位で比べる。生成は、行の総数 82、CI 70、全件だけ 12、終了状態の内訳を、登録した定数と照合する。

### 5.1 定義と数

- **行**は「新しい入口」（入口、入力、設定の組が他の行と違い、自分の期待する終了状態を持つ実行）を1つ持つ。**行内の検査**（KE の対照、基準の子プロセス、再生）と**再利用**（`BE1.result` のように書く）は行に数えない。CI の行は全件だけの行を再利用しない。skip、比較 0 件、期待と違う例外は、合格にも検出にも数えない。行の一意性は、入口ごとの最終の引数で作った署名で見る（別名を値に、欠陥の操作を軸・記録・JSON に適用し、核は名前と窓に、既定値を展開する。数値は最終の binary64 にそろえて直列化する。同じ値を明示する操作は、文字列でも数値でも署名を変えず、CLI のキーの欠落と明示の `null` は別の署名になる。生成が自己検査する）。
- **生成が有理数で証明するもの（各行の `geometry`）:** 換算、単射と最初の衝突の組・出所、\(S_B\)、\(S_K\)、\(S\)、各ビンの \(\lceil U_i\rceil_{64}\)、`complete` か届くかの別（`bound_class`）、食い違いの件数と先頭の例、警告の有無、許容のないときの状態、丸めた KE の問題のゼロ損失（KE の幅で配った値）。登録した定数（2つの台の次数、状態）と照合する。
- **生成が証明できないもの（`expect` の `basis` に出所を書く）:** 有限の上界と許容の大小（BT2a・BT2b は M4b-1 の T2a・T2b の像、BS1t・BS1x は手計算の見積り）、公開の実行の収束（KE の像を現行のコードで実行して、終了状態だけを見た）。
- **キーの台帳:** 生成は、行ごとに具体の `(行, キー, 比較の種類, 対照)` に展開し、重複を拒む。表には、種類と対照ごとの件数と先頭 3 件の例を書く（`check_plan`）。変異の検出先は、その行の台帳にある具体のキーでなければ拒む。CLI の行は、必須のオプション、オプションごとの値の数、登録した値（`k_total`、20、模型の ID、窓、次数）との一致を検査する（表の `cli`）。

| 群 | 行 | 入口 | 行数 | CI | 分母 |
| --- | --- | --- | ---: | ---: | --- |
| BE | BE1、BE1r、BE2、BE3、BE3r、BE4、BE5a、BE5–BE9、BR1、BZ0 | 解析（BE5–BE7、BZ0 は公開） | 14 | 12 | ビット一致と厳密な換算 |
| BT | BT1a–c、BDS、BT2a–c | 解析 | 7 | 7 | 状態、ビット一致 |
| BS | BS1、BS1t、BS1x、BS1d、BS1dl、BS2、BS3、BS3n2、BS4–BS6 | 解析 | 11 | 11 | 両側の規則 |
| BW | BW1–BW4 | 解析 | 4 | 4 | 警告、ビット一致 |
| BN | BN01、BN04–BN25 | 解析（BN16–BN22 は公開） | 23 | 23 | 状態 |
| BC | BC01、BC01s、BC02–BC12 | CLI | 13 | 10 | 再生と状態 |
| BASE | BASE-ms-s、-ms、-ms-s2、-areas、-analytic、-single、-area1、-n13、-cli、-cli-single | 公開、解析、CLI | 10 | 3 | ビット一致 |
| 計 | | | **82** | **70** | |

- 期待する終了状態の内訳: 成功 40、CLI の完了 8、`error` 32、`not_converged` 2（CI だけでは 33、3、32、2）。
- **台帳の件数（表の `counts`。2つの分母を分ける）:**
  - 固定の形の値（スカラーと配列の要素）25,479: ビット一致 `bit_identity` 23,243、打ち切り `truncation` 868、厳密な換算 `exact_conversion` 599、ビンの所属 `bin_membership` 347、終了状態と例外の属性 `state` 268、診断 `diagnostic` 144、上界の挟み込み `bound_bracket` 9、不等式 1。
  - 再帰して比べる記録（大きさが入力で変わる。1組で1件）`record_identity` 469: `single_loss` の全欄、次数ごとの `mesh`、`diagnostics`、警告の列、基準のコードとの全欄の比較、再生、CSV の全セルなど。
  - 違うキーが同じ参照の量を読むことがあるので、独立な条件の数ではない。

### 5.2 行の内容

**BE（成功。どの行も、行内で物理量の対照（§7.1）と比べる）**

| ID | 軸、一次線、端、核、次数 | 期待 |
| --- | --- | --- |
| BE1 | AX_PO、BD_2、BE_E1、K_A1s、\(N=3\) | 成功、`truncated`、\(S=10\)。増加の入力、不等間隔、両側の外が非ゼロ（`outside.low_binding` に 60.5 eV の線のゼロ損失 1/2 と1次、`outside.high_binding` に 72.5 eV の線の高次）。72.5 eV の δ は内部の境の上で、ゼロ損失は `bin[0]` |
| BE1r | BE1 の一次線と端を逆の並びで | 成功。§7.1 の物理量のキーが、`BE1.result` の **BE の表示のビンの添字だけ**を逆にしたもの（`bin[i]` と `bin[6-i]`）とビット単位で同じ。`kinetic_output_edges_eV` と `kinetic_primary` は、BE1 と同じ KE の昇順で、同じ順のまま一致する。外の名前は動かない |
| BE2 | AX_PA、ほかは BE1 | 成功。\(C=1482.25\)、含意されるずれが参照と一致。登録した 90 値（3次数 \(\times\) 9領域 \(\times\) 3成分と、`spectrum` の 9 領域）が `BE1.result` とビット単位で同じ（この組に限る。共通のずれの一般の性質としては登録しない） |
| BE3 | AX_PA、BB_neg、BE_Eneg、K_A3、\(N=2\) | 成功、`truncated`、\(S=15\)。減少の順のビンの一次線、負の BE、不等間隔。一次線の −0.75 eV より低い部分は `outside.low_binding` のゼロ損失。`hole` の負のビンを切らない |
| BE3r | BE3 の一次線を増加の順で（BB_neg_inc） | 成功。§7.1 の物理量のキー、打ち切りのリスト、KE の座標が `BE3.result` と同じ（入力の BE のビン端の配列は、意図して違う）。`input_primary` は入力の順を保つ |
| BE4 | AX_PO、BD_ends、BE_E4、K_A1s、\(N=2\) | 成功、\(S=3\)。最低の端・内部の境・最高の端の上の δ のゼロ損失が `bin[0]`、`bin[1]`、`bin[2]` に 1/4、1/2、1/4。両側の外のゼロ損失は 0 |
| BE5a | AX_H、BD_731、BE_E5a、K_A1s、\(N=2\) | 成功、\(S=6\)。十進の BE で換算が丸まる（`rounding_max_eV` \(>0\)）。δ と端の 73.1 は同じ binary64（合併する）。両側の幅が丸めの分だけ違うので、上界は §7.2 の挟み込みで見る |
| BE5 | 公開。`synthetic_wide`、\(E_0=5926.9\) eV、AX_H、BD_731、BE_E5（減少）、\(N=2\)（全件だけ） | 成功、\(S=40\)。合成の \(h\nu=6000\) eV。`E0_above_5keV` の警告。`comparison` は `not_matched` |
| BE6 | 公開。Al、\(E_0=1000\) eV、AX_D、BD_249、BE_B120、\(N=3\)、`epsrel` \(10^{-7}\)（全件だけ） | 成功、\(S=119\)。実材料。KE の像は M4b-1 の RB1 の入力 |
| BE7 | 公開。`synthetic_wide`、\(E_0=1000\) eV、AX_D、BD_249、BE_EM、\(N=3\)、\(B_{\rm ref}=250.25\) | 成功、\(S=40\)。\(C-B_{\rm ref}=999\) なので `electron_energy_minus_reference_eV` は厳密に 1。核は与えた \(E_0=1000\) eV で呼ばれる（`_call_kernel` を包んでエネルギーの実引数を読む） |
| BE8 | AX_R、BD_723、BE_E8、K_A1s、\(N=2\) | 成功、\(S=6\)。\(C\) が binary64 でない。\(C\) を先に丸めると、6つの KE の端のすべてが変わる（生成が確かめた） |
| BE9 | AX_PO、BD_dup、BE_E4、K_A1s、\(N=2\) | 成功。同じ位置の2本の δ と、端と共通の座標を保つ。ゼロ損失 1/2、1/2、0 |
| BR1 | AX_8192、BB_u、BE_u、K_A1s、\(N=1\) | 成功、`complete_on_display`、\(S=0\)。中の端 \(3u/4\) が KE で \(8192-u\) に丸まる。**ゼロ損失の期待は丸めた KE の問題の値 1/2、1/2**（元の BE の重なり 3/8、5/8 は参照が診断として記録するだけ） |
| BZ0 | 公開。`synthetic_wide`、AX_D、BD_249zero、BE_EM、`"display_support"` | 成功、`complete_on_display`、\(S_B=S_K=0\)、`computed` 0、`orders` は空。`region_loss_upper_eV` は全部 `None`。対照も `"display_support"` |

**BT（打ち切り。両側が一致する入力。解析）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| BT1a、BT1b、BT1c | AX_S、BD_645、K_T1、\(N=2\)、BE_T1_eq／_dn／_up | BT1a: `truncated`、\(S=3\)（高い BE の端 \(-\) 一次線 \(=3a\)。等号は届く）。BT1b: `complete_on_display`、2、上界は全ビン `(0.0, "complete", None)`。BT1c: `truncated`、3。\(C=128\) なので 1 ulp の端が厳密に写る |
| BDS | BT1a の入力で `"display_support"` | 成功、`complete_on_display`、`computed` 3 |
| BT2c | AX_T2、BD_645、BE_T2、K_T2、\(N=4\) | 成功、`truncated`、\(S=50\)。KE の像は M4b-1 の T2c。`bin[2..5]` は有限の上界と \(\tau>0\)、最高の BE のビン `bin[5]` が最大。`bin[0]`、`bin[1]` は `(0.0, "complete", None)`。3つのリストを別々に照合する |
| BT2a | 同じ入力で許容 0.05 | 成功、`bounded_on_display`（M4b-1 の T2a の像）。全ビンの3つの欄（上界、理由、\(\tau\)）が `BT2c.result` とビット単位で同じ。`required_order_estimate` は、同じ許容で実行した KE の対照と同じ |
| BT2b | 同じ入力で許容 \(2.5\times10^{-3}\) | `truncation_bound_exceeds_tolerance`、`stage` は `truncation`、`required_order_estimate` \(\ge5\)（M4b-1 の T2b の像）で、同じ許容の KE の対照の例外と同じ値（R-C の \(S_B=S_K=50\) が、この比較の前提）。例外のビンごとの3つのリストは入力の BE の順で、`BT2c.result` と同じ |

**BS（両側が食い違う入力。設計 §5 第6版。解析。AX_8192、BD_0、K_T1e（\(a=1+2^{-45}\)）、\(N=1\) を基本にする）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| BS1 | 端 BE_c5（設計 §5 の反例） | 成功、\((S_B,S_K,S,N)=(2,1,2,1)\)、`truncated`。唯一のビンは届く側で、上界は有限で正、理由 `None`、有限の \(\tau\)（`(0.0, "complete", None)` でない）。警告 `axis_conversion_support_disagreement`。丸めた KE の実行は完全で上界 0 なので、打ち切りの欄は KE の対照から写せない |
| BS4 | BS1 の軸だけ AX_16（全部の端が厳密に写る） | 成功、\((2,2)\)、`truncated`、警告なし。打ち切りの欄も KE の対照と一致。**`bin[0]` の上界・理由・\(\tau\) が `BS1.result` とビット単位で同じ**（同じ幅 \(2+2^{-43}\)） |
| BS1t、BS1x | BS1 に許容 \(10^{-2}\)／\(10^{-6}\) | BS1t: 成功、`bounded_on_display`、`required_order_estimate` 1、`bin[0]` の上界・理由・\(\tau\) は `BS4.result` とビット単位で同じ。BS1x: `truncation_bound_exceeds_tolerance`、`required_order_estimate` 2（R-C の \(S=2\)。次数 2 で完全）、例外の3つのリストは `BS4.result` と同じ。上界の手計算の見積りは約 \(10^{-4}\)（Codex が現行の KE のコードで \(1.0095\times10^{-4}\) を得た）で、2つの許容はその約 99 倍と約 1/101 |
| BS1d | BS1 で `"display_support"` | 成功、`complete_on_display`、`computed` 2。物理量の対照には整数 2 を渡す |
| BS1dl | BS1d に `order_limit=1` | `order_limit_exceeded`（`requested_order` 2、`order_limit` 1）、評価 0 回。`requested_order` を R-C の \(\max(S_B,S_K)\) と照合する（例外から2つの側の次数は読まない） |
| BS2 | 端 BE_c5rev、核 K_T1（\(a=1\)）（逆向きの反例） | 成功、\((S_B,S_K)=(1,2)\)、`truncated`、上界は有限で正、警告あり |
| BS3 | 端 BE_c5b | 成功、\((8,8)\)、`truncated`。全体の状態は両側で同じで、`bin[0]` だけが違う（端 \(2+2^{-43}\) は \(2a\) に対して BE が大、KE が小）。`bin[0]` は届く側で、上界・理由・\(\tau\) は `BS4.result` の `bin[0]` とビット単位で同じ。`next_reachability_mismatch_bins` は `[0]` |
| BS3n2 | BS3 で \(N=2\) | 成功、`truncated`。同じ食い違いは計算した次数の診断に移る（`computed_relation_mismatch_count` 1。次数 2、`lower`）。`bin[0]` は両側で `complete` |
| BS5 | 端 BE_c5q（\(2.25+2^{-43}\)） | 成功、\((2,2)\)、`truncated`、**警告なし**（次数、状態、3通りの比較はすべて一致）。それでも幅は BE の側が大きく、大きい方を使う |
| BS6 | BS5 の軸だけ AX_16 | 成功。`bin[0]` の上界・理由・\(\tau\) が `BS5.result` とビット単位で同じ。BS5 の KE の対照（幅 2.25）の上界より大きい |

**BW（記録。BE4 の入力で記録だけを変える。§7.1 の物理量のキー、打ち切りのリスト、KE の座標は `BE4.result` と同じ。記録の欄は意図して違う）**

| ID | 記録 | 期待 |
| --- | --- | --- |
| BW1 | REC_NR | 成功。`energy_axis_field_not_recorded`。`not_recorded` は6つの欄（親の1件ずつ）を §3 の順で |
| BW2 | REC_TILT | 成功。`recorded_emission_angle_differs_from_theta`。補正しない |
| BW3 | REC_CORR | 成功。校正のずれと帯電の補正を記録しても、軸はずれない |
| BW4 | REC_NEST | 成功。`not_recorded` は `("source.monochromatic",)` だけ |

**BN（拒否。どの行も核の評価 0 回。書いていない入力は AX_PO、REC_FULL、BD_ends、BE_E4、K_A1s、\(N=2\)。欠陥は表の `defects` に、パスと操作と値で書いた）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| BN01 | AX_8192、BD_eps、BE_c7、\(N=1\)（設計 §7 の反例） | `axis_conversion_not_injective`。`binding_eV` は \((1,\,1+2^{-45})\)、`kinetic_eV` は \((8191,\,8191)\)、`roles` は `("output_edges_eV[1]", "primary.binding_energies_eV[0]")` |
| BN23 | BN01 で δ の強度を 0（BD_eps0） | 同じ（強度 0 の線も座標の集合に入る） |
| BN04–BN06 | `reference` に両方の形の欄／選んだ形の値が `None`／`kind` が `"vacuum_level"` | `invalid_energy_axis`。`reference`・`exclusive`／`reference.kinetic_energy_offset_eV`・`missing`／`reference.kind`・`unknown_kind` |
| BN07、BN08 | `photon_energy_eV` が `"not_recorded"`／`None` で線源の名前が Al Kα（REC_NAMED） | `invalid_energy_axis`。`photon_energy_eV`・`not_recorded_not_allowed`／`missing` |
| BN09 | AX_PA の `kinetic_anchor_eV` が `"not_recorded"` | `invalid_energy_axis`（`reference.kinetic_anchor_eV`・`not_recorded_not_allowed`） |
| BN10 | `energy_axis` と KE の `DeltaPrimary` | `invalid_energy_axis`（`primary`・`requires_binding_primary`） |
| BN11 | AX_S、BD_645、BE_neg_ke | `invalid_energy_axis`（`output_edges_eV`・`negative_kinetic_energy`） |
| BN12–BN14 | `energy_axis` が素の辞書／`photon_energy_eV` が `nan`／`source` が整数 5 | `invalid_energy_axis`（`energy_axis`・`type`／`photon_energy_eV`・`not_finite`／`source`・`type`） |
| BN24、BN25 | `photon_energy_eV` が −1／`core_level` が空白だけ | `invalid_energy_axis`（`photon_energy_eV`・`out_of_range`／`core_level`・`empty`） |
| BN15 | 端が BE_bad | M4a の `SingleLossInputError`、`invalid_output_edges` |
| BN16 | 公開。`energy_axis=None` で BE の一次線 | M4a の `SingleLossInputError`、`invalid_primary`（KE の経路は変えない） |
| BN17 | 公開。`kind` が未知で、窓も \([0,50]\) | `invalid_energy_axis`（軸の検査は M4a の検査より前） |
| BN18 | 公開。BN01 の入力で `model="b1"` | `axis_conversion_not_injective`（単射は模型の検査より前） |
| BN19 | 公開。BS1 の入力（窓 \([1+2^{-45},50]\)、\(E_0=8192\) eV）で `model="b1"` | `invalid_model` |
| BN20 | 公開。同じ入力で \(N=4\)、`max_convolution_points` 4500 | `work_lower_bound_exceeded`（`resource` は `max_convolution_points`、`lower_bound` 4743、`limit` 4500）。食い違いは拒否にならない |
| BN21、BN22 | 公開。同じ入力／端を BE_c5b にして、`max_order` 65 | `order_limit_exceeded`（65、64） |

**BC（CLI。`synthetic_wide`、\(E_0=1000\) eV、\(L_w=20\) Å、`k_total`、`--model first_loss_13a_chain_v0`。共通の引数は `--loss-window-eV 1 50 --weight-length-A 20 --epsilon-momentum k_total` の7要素。書いていなければ `--max-order 3`）**

| ID | 入力 | 期待 |
| --- | --- | --- |
| BC01s | 新スキーマ。AX_D、REC_FULL、BD_249、BE_EMs、`--loss-window-eV 1 8`、`--max-order 2`（CI 用の小さい形） | 完了、終了 0、`truncated`、\(S=5\)。行内: §4.6 の再生、`run.json` の BE の欄、CSV の見出し・行の順・4つの端の列（2次の列を含む）、ゼロ損失の所属、KE の端と \(C\) が参照と一致。累積の損失の警告は出ない（\(16/1000\)） |
| BC01 | 同じ内容で BE_EM、窓 \([1,50]\)、`--max-order 3`（全件だけ） | 完了、\(S=40\)。BC01s と同じ検査。`cumulative_loss_over_E0_gt_0p1` の警告 |
| BC02 | 新スキーマ。AX_DA、REC_NR、BD_249、BE_EMr、`--max-order 1`（全件だけ） | 完了。含意されるずれ 4.25、`not_recorded` の一覧と警告が `run.json` にある。行の順は `high_binding`、`0`…`6`、`low_binding` |
| BC07 | 新スキーマ。設計 §5 の反例（AX_8192、BD_0、BE_c5、`--loss-window-eV 1.0000000000000284 50`、`--max-order 1`。\(E_0=1000\) eV は CONFIG のまま） | **完了、終了 0**、生成物あり。\((S_B,S_K)=(2,1)\)、`truncated`、`bin[0]` の上界は有限で正、警告 `axis_conversion_support_disagreement`。`bin[0]` の3つの欄（有限で正、理由 `None`、有限の \(\tau\)）と警告を、`run.json` と API の再生の両方で読む（核が違うので、BS4 との数値の比較はしない） |
| BC12 | 新スキーマ。AX_D、BD_m1、端は短い形 BE_grid、`--max-order 1`（全件だけ） | 完了、\(S=3\)。展開した端が登録した binary64 と一致。負の `min` を受け付ける |
| BC03 | `measurement.photon_energy_eV` のキーを消す（記録は REC_NAMED） | `error`、`invalid_schema`、終了 2、出力先なし。属性 `detail` は空でない ASCII の文字列（BC05、BC09 も同じ） |
| BC04 | `measurement.photon_energy_eV` が `"not_recorded"` | `error`、`invalid_energy_axis`（`photon_energy_eV`・`not_recorded_not_allowed`） |
| BC05 | 最上位に未知のキー `observed_counts` | `error`、`invalid_schema` |
| BC06 | 設計 §7 の反例（`--max-order 1`） | `error`、`axis_conversion_not_injective`。標準エラーの `binding_eV`、`kinetic_eV`、`roles`。評価 0 回 |
| BC08 | 新スキーマを `xps-single-loss` に（共通の7要素の引数） | `error`、`SingleLossInputError`、`invalid_schema`（引数は完全で、拒否は JSON を読んだ結果） |
| BC09 | 新スキーマで `axis` が `"kinetic_energy_eV"` | `error`、`invalid_schema` |
| BC10 | `measurement.reference` に `binding_anchor_eV: null` を足す | `error`、`invalid_energy_axis`（`reference`・`exclusive`。構造の誤りではない） |
| BC11 | `xpsloss capabilities` | 終了 0。§4.6 の項目。基準の `capabilities` の全体に、許した削除 1 と追加 2 だけを施したものと等しい（履歴がなければ、この比較だけ `not_run`。今の内容の検査は行う） |

**BASE（`energy_axis=None` と既存の入口。基準のコードの子プロセスとのビット一致）**

| ID | 入口と入力 | 期待 |
| --- | --- | --- |
| BASE-ms-s | `multiple_loss_spectrum`、`synthetic_wide`、D1000、E_Ms、窓 \([1,8]\)、\(N=2\)（CI 用の小さい形） | 成功。BASE-ms と同じ検査 |
| BASE-ms | `multiple_loss_spectrum`、`synthetic_wide`、D1000、E_M、\(N=3\)（全件だけ） | 成功。全欄が基準とビット単位で同じ。基準にない欄はスペクトルの結果の直下の `energy_axis` だけで、値は `None`。行内: `energy_axis=None` を明示した呼び出しも同じバイト列 |
| BASE-ms-s2 | 同じで `multiplier="inverse_window_area"`（全件だけ） | 同じ |
| BASE-areas | `multiple_loss_order_areas`、\(N=3\)（全件だけ） | 同じ。行内: `energy_axis` を渡すと `TypeError` |
| BASE-analytic | 解析の入口、K_A3、D1000、E_A3、\(N=3\) | 同じ（`panels_by_order` も比べる） |
| BASE-single、BASE-area1 | `single_loss_spectrum`／`loss_window_area`（全件だけ） | 同じ |
| BASE-n13 | `multiple_loss_spectrum`、\(N=4\)、`max_convolution_points` 4500（M4b-1 の N13） | `work_lower_bound_exceeded`、評価 0 回。型、`reason`、`vars(exc)` の全部の旧い属性、文言が基準と同じ。明示の `energy_axis=None` でも同じ例外（別のキーで照合する） |
| BASE-cli | `xps-multiple-loss`、旧スキーマ、D1000、E_M（全件だけ） | 完了（両側とも終了 0）。2つの CSV はバイト単位、`run.json` は時刻・所要時間・パス・版を別に確かめて、残りを全部比べる（`energy_axis` のキーがどこにもない） |
| BASE-cli-single | `xps-single-loss`、旧スキーマ、共通の7要素の引数（全件だけ） | 同じ |

- **基準のコード:** main のコミット `8c33d2939dab866decc37fdf9663f5ac2e29b78d`（M4b-1 の完了。3つの OS の CI が通った `e7fa6a3` と `src/` が同じ）。M4b-1 は未公開でタグがないので、完全な SHA で固定する。`energy_axis` を足す対象がこのコミットにあるので、M4a・M4a-2 の基準（`v0.3.0`）では足りない。
- **方法**（M4a-2 の V1、M4b-1 の M4A と同じ）: 1つの子プロセスで `git archive <完全な SHA> src/xpsloss`（同梱のデータを含む）を使い、選んだ BASE の行をまとめて動かす。親と同じ実行ファイルと依存の版、子のパッケージと実装のモジュールの実際のパスを確かめる。チェックアウト（作業ツリーの `.git` ファイルを含む）で、コミットがない・一致しない・子の読み込み先が違えば**失敗**。**展開した sdist で `.git` がなければ、BASE の比較は `not_run/sdist_without_git_history` と記録し、合格に数えない**（M4b-1 の `tests/m4b1_states.py` と M4a-2 の扱いに合わせる）。
- **比べ方:** 構造、型、辞書のキー、警告の順、配列の dtype・形・バイト列、スカラーのビット。BASE で許す追加は、スペクトルの結果の直下の `energy_axis is None` だけ。版の欄は、各プロセスが自分の `__version__` を報告することを照合する。
- 既存の試験と参照（M4a の 100 行、M4a-2 の行、M4b-1 の 109 行）は、§9.1 の変更のほかは変えずに合格すること。

### 5.3 CI と全件

- **CI（70 行）**と**全件だけ（12 行）:** BE5、BE6、BC01、BC02、BC12、BASE-ms、BASE-ms-s2、BASE-areas、BASE-single、BASE-area1、BASE-cli、BASE-cli-single。BC01 と BASE-ms は、同じ入口・材料・換算の小さい形（BC01s、BASE-ms-s）を CI に置いた。全件だけの行も、全件と変異の対象である。CI の行は全件だけの行を再利用しない（生成が検査する）。
- CLI の行は同じプロセスの中で `main()` を呼ぶ。基準のコードの側は1つの子プロセスで行う。

## 6. 参照値（実装担当とは別の担当。R-C）

参照の担当は、実装のコードと結果を見ずに作る。実装担当は参照値を見ない。参照は凍結して SHA-256 を付録に記録してから、実装と比べる。

- **生成:** `benchmarks/numerics/m4b2/generate_references_m4b2.py`（標準ライブラリの `fractions` だけ。`xpsloss` を読み込まない）が `references_m4b2.json`（ASCII、LF）を書く。`--check` でバイト単位に再生成する。入力は凍結した表（SHA-256 を照合してから読む）。
- **対象の行:** BE の 14 行、BT の成功の 6 行、BS の成功の 9 行、BW の 4 行、BC01s・BC01・BC02・BC07・BC12（換算、所属、両側の台）、BN01・BN23・BN18・BC06（衝突の組と出所）、BT2b・BS1x・BS1dl（両側の台の次数）。
- **中身:** 入力の座標ごとに1回だけ \({\rm RN}(C-B)\) を作り（16進）、その座標を使う出力（KE の端の各要素、KE の一次線の各要素）へ配る。\(C\) の既約分数と \({\rm RN}(C)\)、`rounding_max_eV`、含意されるずれ、`primary_reference_kinetic_eV`、`electron_energy_minus_reference_eV`。物理量の対照の入力（KE の昇順の端と一次線、渡す `max_order`）。§2.3 の \(S_B\)、\(S_K\)、\(S\)、各ビンの \(\lceil U_i\rceil_{64}\)、`complete` の別、食い違いの件数と先頭の例、警告の有無、KE の側との幅の差の最大。領域ごとのゼロ損失（**丸めた KE の幅で配った値**の厳密な有理数。ビンの一次線で丸めが起きる行では、元の BE の重なりを診断として別に記録する）。BE のビンと KE のビンの対応。2つの端のハッシュ（BC の行）。
- **独立性:** 製品の変換の関数を呼ばない（往復で確かめない）。丸めは `float(Fraction)`、所属は有理数の比較で決める。表の生成の中の計算（§5.1）とは別の担当の別のコードで、2つが食い違えば、実装と比べる前に止めて原因を記録する。
- **自己検査と既知の誤答:** AX_D と AX_DA の KE が一致すること。設計 §5・§7 の反例。BE8 で \({\rm RN}({\rm RN}(C)-B)\ne{\rm RN}(C-B)\)（6つの端）。X25・X36・X41 の誤った \(C\) での KE の端。BS1 で KE の側だけの判定（\(S=1\)、完全）。BR1 の 3/8 と 1/2。
- **参照にしないもの:** スペクトルの値（次数ごとのビン）と上界の値の独立な数値の参照は置かない。値は物理量の対照とのビット一致で、上界は §7.2 の規則で確かめる。KE の値は M4b-1 の検証（R-A、R-B）に依る。M4b-1 の参照値をここで読み直さない。
- **「値なし」:** 参照のキーはすべて厳密な有理数か binary64 で、\(u_r=0\)。`no_reference` は 0 件のはずで、1つでも欠ければ、その行は「参照不足」と記録して合格に数えない。

## 7. 比較のキーと許容幅

### 7.1 キー、写像、2つの対照

- **キー:** `order[n].<領域>.<total|no_hole|hole>`、`zero_loss.<領域>`、`spectrum.<領域>`（領域は `outside.low_binding`、`bin[i]`、`outside.high_binding`）、`order[n].area.<列>`、`truncation.*`、`energy_axis.<欄>`。読み場所は §4.3。
- **写像:** 増加の入力なら BE の `bin[i]` は KE の `bin[n_b-1-i]`、減少の入力なら KE の `bin[i]`。`outside.high_binding` は KE の `outside.low`、`outside.low_binding` は KE の `outside.high`。
- **物理量の対照（KE の対照とのビット一致）:** 同じ入口を、R-C の KE の端と一次線、同じ \(E_0\)・\(L_w\)・\(m\)・設定で、**BE の実行が実際に計算した次数**で実行する（`"display_support"` で次数が違いうるので、KE には BE の側が決めた正の整数を渡す。両側とも 0 のときだけ KE も `"display_support"`）。許容は外す。比べる対象は、各次数の3成分の値・返す誤差・4つの誤差の内訳（\(n_b+2\) 領域）、面積・面積の誤差・感度・感度の誤差・`area_zero_kind`（11 列）、面積の誤差の内訳（5 \(\times\) 11）、`zero_loss`・`spectrum` と誤差、\(\hat r\) と区間、`convention_ratio`、記録として `single_loss` の全欄、次数ごとの `mesh`、`neglected_total_area`、`diagnostics`（軸に関わらない部分）、警告（BE の識別子と `truncated` を除いた識別子の列と順）。
- **打ち切りの対照（両側の規則。R-C）:** `truncation.support_order`、`state`（**許容のない行。許容を与えた成功の行 BT2a・BS1t の `state` は、上界の数値に依るので R-C でなく `expect` の `bounded_on_display` と照合する**）、`computed`、`requested`、ビンごとの `complete` か届くかの別、`energy_axis` の §2.3 の欄、警告 `axis_conversion_support_disagreement` の有無は、R-C の値と照合する。**状態に依る欄（`truncated` の警告、`requested`、許容、必要な次数、上界・理由・最適の \(\tau\)）を、物理量の対照との無条件の一致に入れない。**
  - **`complete` か届くかの別は、実際に返った3つの欄から判定する。**`complete` は3つの欄がちょうど `(0.0, "complete", None)`、届くビンはそれ以外。
  - **ビンごとに決める。**両側が届き、上へ丸めた幅が KE の側と同じビン（表の `width_equals_kinetic`）は、上界・理由・\(\tau\) を KE の対照とビット単位で照合する。
  - 両側が届き、幅が広がるビン（BE5a、BE5、BE8、BS5 など）は、上界を §7.2 の挟み込みで照合する。
  - KE の側では届かないビン（BS1・BS1t・BS3 の `bin[0]`、BC07）は、3つの欄が「有限で正、理由 `None`、有限の \(\tau\)」であることを照合する（`(None, "no_usable_tau", None)` では通らない）。**数値のビット一致は、BS1・BS1t と BS3 の `bin[0]` を、端が厳密に写る対照の行 BS4 と比べるものに限る**（同じ核・一次線・次数で、1次のパネルも同じ）。BC07 は核が違う（公開の `synthetic_wide`）ので、BS4 と比べない。
  - **状態と警告:** `truncated` の警告は、`truncation.state == "truncated"` と一対一で照合する。許容は入力と照合し（許容を外した KE の対照と `diagnostics.effective_settings` を再帰で比べるとき、許容の欄だけは入力との照合に戻す）、許容のもとでの終了状態は `expect`、幾何は R-C、数値は対照で見る。
  - **失敗の行で使う R-C:** BT2b は \(S_B=S_K=50\)（同じ許容の KE の対照の例外と比べる前提）、BS1x は \(S=2\)（`required_order_estimate`）、BS1dl は \(\max(S_B,S_K)=2\)（`requested_order`）。
- 合否は「行、キー、比較の種類、対照」の単位で数える。登録した参照のキーの未使用と重複の使用は失敗にする。

### 7.2 許容幅

| 種類 | 幅 |
| --- | --- |
| `exact_conversion`（\({\rm RN}(C)\)、KE の座標、含意されるずれ、\(B_{\rm ref}\) の KE、\(E_0\) との差、`rounding_max_eV`） | binary64 のビット一致（符号つきのゼロを含む）。`conversion_constant_exact` は有理数として等しいこと |
| `bin_membership`（ゼロ損失） | 参照が 0 の領域は値も誤差も厳密に 0。ほかは \(|y-r|\le4\epsilon|r|\)。**\(r\) は丸めた KE の幅で配った厳密な有理数** |
| `bit_identity`、`record_identity`（KE の対照、別の行、基準のコード、再生） | ビット単位（配列のバイト列。符号つきのゼロを含む） |
| `truncation` | 整数、状態、`complete` の別、件数、添字は完全一致。`region_loss_upper_eV` は binary64 のビット一致 |
| `bound_bracket`（両側が届き、幅が広がるビンの上界） | \(b_K\le b\le b^{\uparrow}\)。\(b_K\) は KE の対照の上界。\(b^{\uparrow}\) は、**同じ1次のパネルに対し、KE の対照の最適の \(\tau_K\) を固定し、R-C の幅 \(\lceil U_i\rceil_{64}\) で、M4b-1 と同じ順に外向きに評価した1つの候補**（読む場所は §4.5）。その候補が新しい幅で表せなければ、有限のほかの候補の最小を上限にし、どれもなければ理由を照合する。厳密なゼロ、正の下へのあふれ、`no_usable_tau` は M4b-1 の分岐のまま |
| 状態、理由、属性、警告の識別子、記録の欄 | 完全一致。属性の数値は binary64 のビット一致 |

**上界の上限に、幅の差の式を使わない理由（Codex の反例）:** 第2稿の上限 \(b_K(1+2\tau_{\max}\Delta U)+4\epsilon b_K\) は、binary64 の一般の上限にならない。K_T1e、\(m=2\)、\(N=1\)、BE 0 の δ の強度を \(e^{-600}\)、表示の端を 0 と `0x1.0000000000085p+1` eV にすると、\(C=8\) では KE の幅が `0x1.0000000000084p+1`、\(C=4\) では元の幅が厳密に写り、\(\Delta U=2^{-51}\) eV。同じ \(\tau\) で \(b_K=2.6756594734474793\times10^{-265}\)、\(b=2.6756594734477834\times10^{-265}\) で、相対の増分 約 \(1.14\times10^{-13}\) は、旧い式の許す \(3.64\times10^{-14}\) を超える（対数の和の 1 ulp は、指数に戻すと相対 \(4\epsilon\) より大きくなりうる）。新しい上限は、実数では \(\min_jf_j(w)\le f_{j_K}(w)\) から成り立ち（最適の \(\tau\) が替わっても通る）、同じ演算の順で評価するので丸めにも依らない。下限は、幅についての単調性から成り立つ。登録した 9 キーでは、Codex の試算で最適の候補の交替はなかった。

**\(4\epsilon\) を「元の BE の厳密な値」との比較に一律に使わない理由（Codex の例）:**

- \(C=8192\)、一様な一次線の BE のビン \([0,2u]\)（強度 1）、表示の端 \([0,3u/4,2u]\)。合併は単射だが、中の端は KE で \(8192-u\) に丸まり、低い BE のビンのゼロ損失は元の BE では 3/8、丸めた KE では 1/2 になる（BR1）。
- BE 0 の δ に、\(a=1\)、\(b=1+7u/8\) で両端が 0、中央が頂点、面積 1 の三角形の連続な核を使うと、表示の上端 \(1+3u/4\) までの1次の面積は元の BE で 47/49、丸めた KE では 1。\(nb\) をまたぐ効果も「丸めの水準」ではない。
- したがって期待は丸めた KE の問題の値とし、元の BE の重なりは診断にだけ使う。換算が厳密な行（BE3 など）では、2つは一致する（生成が確かめる）。

### 7.3 判定と数え方

- このキーの大半は、厳密な参照かビット一致なので、**合格**か**不合格**になる。前提が成り立たない行（基準の子プロセスが作れない、参照が欠ける、物理量の対照が未収束）は**未実行**か**参照不足**として別に数え、合格に含めない。必須の行の未実行を合格にしない。
- 返す誤差の推定の検査（M4b-1 の §7.3）は、この段では行わない（独立な数値の参照がない）。
- **別々に数えるもの:** 合格、不合格、`no_reference`、`not_run`（理由つき。`sdist_without_git_history` を含む）、実装の予期しない未収束、上流の失敗、参照の経路の失敗、参照不足の行。未使用・重複・欠測・対照の失敗は、台帳の分母と突き合わせる。

## 8. 時間の上限

- **CI の増分:** 各代表の環境で 60 秒以内（製品の実行、参照の読み込み、対照、基準の子プロセス、再生を含む経過時間）。**下は Codex の実測に基づく見積りで、M4b-2 の実装の実測ではない**（この PC、Windows、Python 3.12.10、直列。BE の入口を丸めた KE の入力で代用し、CLI は3つのファイルを書き、再生はその manifest から API を呼び、基準は完全な SHA の `src/xpsloss` を展開した子1つ）。
  - 第2稿の CI の公開の行は、公開に当たる計算 21 回で 120.50 秒、一式の経過で **128.64 秒**（上限の2倍を超える）。行ごとには BE5 9.49、BE7 13.16、BZ0 10.83、BC01 14.26、BC02 8.40、BC07 7.00、BC12 9.57、BASE-ms 20.20、BASE-areas 10.95、BASE-cli 16.65 秒。
  - この第3稿の CI の公開の行（BE7、BZ0、BC07、BC01s、BASE-ms-s と、評価 0 回の行）は、公開に当たる計算 11 回で、子の起動を含めて **38.99 秒**。解析の行（BE・BT・BS・BW の成功と打ち切りの超過）は対照を含めて約 1.83 秒。**見積りは約 41 秒**で、BE の実装、R-C の読み込み、台帳の全比較（未計測）に約 19 秒を残す。
- **60 秒を超えたときの順（表の `ci_reduction_order`）:** BZ0（10.83 秒。**空の台を最後まで通す、ただ1つの登録した行**なので、回すと空の台の実行は全件だけになる）、BE7（13.16 秒。回すと、X42 と、X25 の3つ目の枝（\(\phi\) を落とす）の検出も全件だけになる）、BC01s の順に、各環境の実際の増分が 60 秒を超えたときだけ1行ずつ全件だけに回し、回すたびに測り直す（CI／全件だけは 70／12、69／13、68／14、67／15）。BC07（反例の CLI）、BASE-ms-s（基準）、BASE-analytic は残す。回した行は全件と変異の対象に残す。それでも超えれば、旧版を残して再登録する。
- **全件:** 30 分以内（見積りは 5 分以内。第2稿の公開の行の 128.64 秒に、BE6 の Al の2回の実行 約 30–40 秒などが加わる）。**1行:** 120 秒を超える行は登録しない。
- **参照の生成:** 5 分以内（有理数の計算だけ）。**全変異:** 必須 14 件の正常版と変異版の対（枝ごと）を 60 分以内。

## 9. 前の凍結した段への影響と、変異

### 9.1 前の凍結した段の表と試験

試験のコードを読んで、M4b-2 の変更が当たる所を挙げた（Codex も読んで確かめた。コードの読解による判定で、実行はしていない）。

| 段と行 | 試験の検査 | M4b-2 の変更 | 扱い |
| --- | --- | --- | --- |
| **M4b-1 C11**（`tests/m4b1_states.py:1624`。`UNSUPPORTED` との比較） | `multiple_scattering.unsupported_physics` が 10 個の列と**順つきで完全一致** | `binding_energy_axis` を外す | **凍結した段の試験の変更（1）。**登録した 1 要素だけを除いた 9 要素との順つきの一致にする（M4a の `LATER_SUPPORTED` と同じ形）。凍結した表は変えない |
| M4b-1 の単体の試験 `tests/test_m4b1_cli.py:135` | `binding_energy_axis` が `unsupported_physics` に**ある**こと | 同上 | **試験の変更（2）。**「ない」ことの検査に変える |
| M4b-1 C09、`test_rejection_records` | 旧スキーマで `axis` が `"binding_energy_eV"` なら `invalid_schema` | 旧スキーマは KE だけのまま | 変えずに済む（§4.6） |
| M4a の `defect` の行、`tests/test_m4a_cli.py` | `xps-single-loss` が BE の `axis` を拒否 | 変えない | 変えずに済む |
| M4b-1 C11 の1・2つ目、M4a V08 の `LATER_BLOCKS`、M4a-2 の `test_m4a2_records.py:678` | `supported` が基準 \(+\) `xps_multiple_loss`。旧いブロックは基準（`v0.3.0`）と完全一致。追加のブロックは比較から外す | `supported` は変えない。変えるのは `multiple_scattering` だけ（旧いタグの側にない） | 変えずに済む |
| M4b-1 C01 | 設定のキーが 13 個、出力のファイルが3つ、旧い CSV の見出し、再生 | 設定は足さない。KE の `run.json`・CSV・失敗の JSON に欄を漏らさない | 変えずに済む |
| M4b-1 の M4A の群 | 内包した旧い結果の全欄を `v0.3.0` と比べる | `SingleLossSpectrumResult` は変えない | 変えずに済む |
| 旧い入口を名前で呼ぶ試験、M4a-2 のキーワード専用の結果の欄の検査 | 旧い型と署名 | 保つ | 変えずに済む（公開の名前の全体や、多重損失の関数の署名の全体を固定する試験は見当たらない） |

**凍結した段の試験の変更は上の2か所だけ**で、利用者の了解の対象もこの2か所に限る（了解は凍結のときに統括担当が取る）。新しい BC11 は、基準の `capabilities` の全体に、許した削除 1（`binding_energy_axis`）と追加 2（`energy_axes`、`measured_spectrum_schema`）だけを施して比べる。`binding_energy_axis` を残す案は、対応した機能を「未対応」と表示し続けるので採らなかった。

### 9.2 変異の表（必須 14）

レビュー担当が、各変異を本番のコードに1つずつ入れる。設計 §10 の X25、X26、X27、X32（設計由来 4）と、この段で足した X34–X43（追加 10）を、すべて必須にした。検出先は表の `mutations` に、枝ごとに `(行, キー, 比較の種類)` と、結果に届かない枝の予告した型・理由で書いた（生成が、台帳にあるキーであることを検査する）。下は実装の前に成り立たせる検査で、検出の実績ではない（未測定）。

| ID | 変異（枝） | 登録した行と、正常版から区別する観測 |
| --- | --- | --- |
| X25 | (1) \(E=B-C\) (2) \(E=C+B\) (3) \(\phi\) を落とす（\(C=h\nu\)）。3つを別々に走らせる | (1) BE1 が `invalid_energy_axis`・`negative_kinetic_energy` になる（予告した型と理由）。(2) BE1 の KE の端、`zero_loss.outside.low_binding`、`order[1].bin[4].total`（§2.2 の狭義の減少を検査する実装では `axis_conversion_not_injective`）。(3) BE7 の `conversion_constant_eV`、KE の端、`electron_energy_minus_reference_eV`（軸の共通のずれはビンの値を変えないので、記録で見る） |
| X26 | (1) 出力を入力の向きに戻さない (2) 外の名前を入れ替えない | (1) BE1 の `zero_loss.bin[0]`、`order[1].bin[4].total`（損失 3–5 eV のビン。逆の並びでは損失が届かないビンに当たる）、BE1r。(2) BE1 の `zero_loss.outside.low_binding`（1/2 のはずが 0）、`order[3].outside.high_binding.total` |
| X27 | 内部の境の上の δ を高い BE 側のビンへ入れる | BE4 の `zero_loss.bin[1]`、`bin[2]`（1/2、1/4 が 0、3/4）、BE1 の `zero_loss.bin[0]` |
| X32 | 丸めた KE の差だけで台を判定する | **両方の版が成功する**ので、終了状態では見えない。BS1 の `support_order`（2 が 1）、`state`（`truncated` が完全）、`bin[0]` の別と3つの欄（有限で正が `(0.0, "complete", None)`）、`support_order_binding` |
| X34 | 2回丸める（\({\rm RN}({\rm RN}(C)-B)\)） | BE8 の6つの KE の端のビット列と `rounding_max_eV` |
| X35 | 合併の検査だけを弱める（配列ごとの単射は残す） | BN01 は成功、BC06 は完了になる（登録した型・理由・組・`roles` が消える。丸めた KE の像は成功する）。ほかの例外は検出に数えない |
| X36 | 対応点の和の符号（\(C=E_\star-B_\star\)） | BE2 の `conversion_constant_eV`（1314.25）と KE の端 |
| X37 | 欠けた \(h\nu\) を線源の名前から 1486.6 eV で補う（スキーマの検査より前に入れる）。(1) API の明示の `None` (2) CLI のキーの欠落 | (1) BN08 (2) BC03。どちらも REC_NAMED で、登録した型・理由が消える |
| X38 | (1) `"not_recorded"` の警告を出さない (2) \(h\nu\) の `"not_recorded"` を有限の値 1486.6 に置き換えて通す | (1) BW1 の警告の列（`not_recorded` の一覧が残っていても落ちる。一覧が空になることは求めない） (2) BN07 が成功する（偶発の `TypeError` は検出に数えない） |
| X39 | 値は並べ直すが、`truncation` のリストを並べ直さない。(1) 上界 (2) 理由 (3) \(\tau\) | BT2c の `bin[5]` の上界、`bin[5]` の理由（`truncation.bound_reason[bin[5]]`。`None` が `complete` になる）、`bin[5]` の \(\tau\) を別々に直接読む |
| X40 | (1) 台の次数は両側で取るが、ビンごとの幅は KE だけ (2) 計算した次数の診断を省く (3) 比較が食い違ったときだけ大きい幅を使う | (1) BS3 の `bin[0]` の別、3つの欄、`region_loss_upper_eV`（全体の状態は両方の版で `truncated`）。(2) BS3n2 の件数と先頭の例。(3) BS6 と BS5 の上界のビット一致 |
| X41 | 対応点の含意されるずれを補正に使う（\(C=h\nu\)） | BE2 の `conversion_constant_eV` と KE の端、BC02 の `implied_offset_eV` と KE の端（記録だけ正しい変異にも効く） |
| X42 | \(E_0\) を \(C-B_{\rm ref}\) から決める | BE7 の `order[1].bin[5].total`（Codex が、KE の対照で \(E_0\) を 999 eV にするとビットが変わることを確かめた）と、核に渡ったエネルギーの実引数 |
| X43 | 全体の端の上の δ を外の領域へ入れる。(1) 低い端 (2) 高い端 | BE4 の `zero_loss.bin[0]`・`outside.low_binding`／`bin[2]`・`outside.high_binding` を別々に |

**検出の数え方**は M4b-1 の §9 と同じ: 登録した入力で正常版が通り、変異が実際に走り、登録したキーか、予告した型・理由の差が出て、初めて1件。正常版の合格と変異版の不合格をそろって必須にする。枝のある変異は枝ごとに結果を残し、全部の枝が落ちたときに検出済みとする。skip、比較 0 件、期待と違う例外、試験器の例外、変異が入っていない場合、参照なし、正常版の失敗、登録の外の検査だけの失敗は、検出に数えない。報告は必須 \(k/14\) とし、設計由来 4 と追加 10 の内訳、枝ごとの成否、生存、未実行、参照不足を別々に示す。

## 10. 凍結の前の確認

**下書きの担当が確かめたこと（2026-11-28）:**

- 設計 §5 と §7 の2つの丸めの反例は、`fractions` と `float` だけのスクリプト（`work/consult/m4b2/prereg/counterexamples_check.py`。パッケージを読み込まない）で成り立つ。表の生成も、BS1（\((2,1)\)）、BS2（\((1,2)\)）、BS3（\((8,8)\)、`bin[0]` だけ届き方が違う）、BS3n2（計算した次数の食い違い 1 件）、BS5（食い違いなしで幅だけ違う）、BS4・BS6（厳密に写り、幅が BS1・BS5 と同じ binary64）、BR1（1/2 と 3/8）、BE8（6つの端）を有理数で確かめる。
- `_fold.Geometry.mass` の境の規則（§2.4）。設計 §7 の文と食い違わない。
- 材料の2つのファイルと M4b-1 の表の SHA-256 は現在のバイト列と一致する。基準の `8c33d29` と `e7fa6a3` の `src/` の差は空。
- 表の生成は、82 行・CI 70・全件だけ 12、終了状態の内訳、最終の引数で作った署名の一意性、語彙、新しい理由ごとに1行以上、外した理由が残っていないこと、再利用の向き、CI を回した各段階の件数、台帳の重複、変異の検出先が台帳にあること、CLI の引数（必須のオプション、値の数、登録した値）を例外で検査する。`--check` はバイト単位で一致する。
- 公開の行の KE の像（BE5、BE7、BZ0、BC12）は、現行のコードで成功する（終了状態と時間だけを見た。`timing_ke_only.py`、`haxpes_row_check.py`、`state_check_draft2.py`）。Codex は、ほかの成功の行の KE の像の状態、BC07 の KE の像の成功、修正した CLI の引数、BE1 と BE2 の 90 値のビット一致、X42 を確かめた。これらは M4b-2 の合格の件数に入れない。

**Codex の1回目の確認（判定「修正後に凍結」）への対応:**

| 指摘 | 対応 |
| --- | --- |
| 1（blocker）BC08 と BASE-cli-single の引数 | 共通の部分を `k_total` までの7要素にした。生成が、全部の CLI の行で必須のオプションと値の数を検査する |
| 2（blocker）(d) と KE の対照の全キーの一致は両立しない | §7.1 を物理量の対照と打ち切りの対照に分けた。§2.3 の規則を次数の計画と上界の両方に入れる契約にした |
| 3（major）丸めたビンの一次線 | 返す数値は丸めた KE の問題の値と明記し（§2.1）、ゼロ損失の期待を KE の幅で配った値にした。BR1 を足し、§7.2 に2つの例を書いた |
| 4（major）読み場所 | §4.3（最上位の `outside[side]` の4欄、読み場所の表、内部の KE の索引のまま残るもの）、§4.6（BE の CSV の作り方） |
| 5（major）`run.json` の形と再生 | §4.6（KE の manifest に `energy_axis` を出さない、`input_primary` の形、再生に使う欄、全セルの照合、2つのハッシュ） |
| 6（major）構造と値の判定の順 | §3 の「判定の順」と値の規則（深い写し、`out_of_range`、`empty`）、§4.2。BN24、BN25、BC10 |
| 7（major）生成器 | 不正な軸は「正しい軸 \(+\) 1つの欠陥の操作」で書き、記録の辞書から説明のキーを除いた。BC03 は REC_NAMED。有理数で証明する幾何と、出所つきの終了状態を分けた。署名は別名と既定値を展開する。Al の実行用の設定のハッシュも照合する |
| 8（major）キーの台帳 | `(行, キー, 比較の種類, 対照)` に展開し、固定の形の値と再帰の記録の分母を分けた（§5.1） |
| 9（minor）片方向だけの仕様 | BE3r、BE9、BN23、BZ0、BW3、BW4、BC12 を主行として足した |
| Q1（台の判定） | 統括担当の決定で (d)。§2.3、BS の群、BN19–BN22、BC07 |
| Q5、Q9、Q10、Q11 | 必須 14。BE1・BE2 の 90 値をビット一致に。sdist は `not_run`。BASE-n13 を足した |

**Codex の2回目の確認（判定「修正後に凍結」）への対応:**

| 指摘 | 対応 |
| --- | --- |
| 1（blocker）BC07 と BS4 は核が違う | §7.1 のビット一致の対象を、BS1・BS1t と BS3 の `bin[0]` と BS4 に限った。BC07 は3つの欄と警告を `run.json` と API の再生で読む |
| 2（major）`bound_bracket` の右辺 | §7.2 を「下限 \(b_K\)、上限は \(\tau_K\) を固定した1候補」に替え、ビンごとに決めることにした。読む場所は §4.5（既存の非公開の `_order_bounds.truncation` と `tail_model`）。反例を書いた |
| 3（major）次の次数の例 | §2.3 の 7。`first_mismatches` に次の次数の例を入れ、両方の軸の座標、入力の添字、比較の種類を持たせ、順を固定した。表の `geometry` も同じ形 |
| 4（major）打ち切りの成功の観測 | §7.1 と台帳: `truncated` の警告と状態の一対一、3つの欄からの判定、BT2a の全ビンと必要な次数、BS1t、BS3 の理由と \(\tau\)、失敗の行の R-C、許容の欄の扱い |
| 5（major）CI の時間 | BC01 と BASE-ms を全件だけにし、小さい形 BC01s・BASE-ms-s を足した。BASE-cli、BC12、BE5、BC02、BASE-areas を全件だけに移した。82 行、CI 70（§8） |
| 6（minor）署名と台帳 | 署名を最終の引数で作る。`check_argv` に値の照合。BE1r・BE3r・BW の範囲。`detail`、終了 0、明示の `None`、評価 0 回のキー |
| 7（minor）変異 | X35、X37、X38、X39 の予告と読み場所を直した |
| 設計 §5 の文案 | 統括担当が設計に入れた。§2.3 の 1・2・6・7 と食い違わない |

**ここで決めた事項（設計が名前や値を書いていない所。第2稿・第3稿で加えたものを含む）:** 型の名前と欄、記録の欄に既定値を置かないこと、軸の型は構築のときに検査しないこと。記録の辞書の形と値の規則（§3）、`charge_correction.applied` が `"not_recorded"` のときの形、`"not_recorded"` の一覧の順。`primary_reference_binding_eV` は数値を必須にすること。単射の集合の範囲と報告する組・出所。負の KE を `invalid_energy_axis` にし、記録の値の負は許すこと。結果の `energy_axis` のキー、最上位の `outside[side]` の4欄、`single_loss` を KE のまま持つこと。例外の属性と `detail` の語彙、警告の4つの識別子と順。新スキーマの欄の名前、`axis` を BE だけにすること、判定の順、短い形の端の展開。BE の CSV の見出しと4つの端の列、`run.json` の BE の欄、KE の実行の出力を変えないこと。`capabilities` の変更。基準のコード。上界の挟み込みの形と、照合器が読む場所（§7.2、§4.5）。`first_mismatches` の例の欄の名前（§2.3）。BS1t・BS1x の許容（\(10^{-2}\)、\(10^{-6}\)）。CI 用の小さい形の入力（窓 \([1,8]\) eV、\(N=2\)、端 E_Ms）と、全件だけに回す順。

**Codex の確認の記録（3回）:**

| 回 | 対象 | 判定 | blocker | 記録 |
| --- | --- | --- | ---: | --- |
| 1 | 第1稿 | 修正後に凍結 | 2 | `work/consult/m4b2/prereg/codex_prereg_review.md`（major 6、minor 1） |
| 2 | 第2稿 | 修正後に凍結 | 1 | `codex_prereg_review2.md`（major 4、minor 2。1回目の指摘はほぼ解消） |
| 3 | 第3稿 | 修正後に凍結 | 0 | `codex_prereg_review3.md`（技術の major 1、minor 3。§4.5 の約束は実装担当として受け入れられる） |

3回目の「残る指摘」への対応（第4稿）: (1) BT2a・BS1t の `truncation.state` の対照を `expect` にした（§7.1。行の数と種類ごとの分母は同じ）。(2) 署名の数値を最終の binary64 にそろえ、数値の同じ値の操作を自己検査に足した。(3) BE1r の「逆」を BE の表示のビンの添字に限ると書いた。(4) §8 に、BZ0 が空の台を最後まで通すただ1つの行であることと、BE7 を回したときの X42・X25 の3つ目の枝を書いた。Codex の試算は現行の KE を対象にしたもので、M4b-2 の実装の後の照合、CI の時間、変異の検出は未測定である。

**統括担当の凍結の記録:** 2026-11-28 に利用者が凍結を了解した。同時に、設計 §5 第6版の BE の台の判定（両側で求めて弱い側の主張を返し、`axis_conversion_changes_support` で拒否しない）と、§9.1 の M4b-1 の試験2か所（`tests/m4b1_states.py` の C11 の `UNSUPPORTED` の完全一致、`tests/test_m4b1_cli.py::test_capabilities`）を、登録した削除 `binding_energy_axis` だけを受け入れる形に変えることを了解した。表は 82 行（CI 70、全件だけ 12）、必須の変異 14。4回目の Codex の確認は省き、統括担当が第4稿の修正を直接確かめた（`--check` を `PYTHONHASHSEED` 11 と 271 で一致、BT2a・BS1t の `truncation.state` の対照が `expect`、§10 の記録）。ハッシュは付録、凍結のコミットはその記録 `work/consult/m4b2/prereg/freeze_record.md`。

## 11. 凍結の後の変更

（凍結の時点では空。）凍結の後に期待する状態、許容幅、行、変異を変えるときは、利用者の了解を得て、次を記録する: 変える行とキー、凍結した期待と変えた期待、理由、何が起きたか、原因、確かめたこと。凍結した表はバイト単位で変えずに残し、試験は表を読んだ後に変えた値だけを置き換える（M4a の `POST_FREEZE_EXPECT` と同じ）。参照値を見た後の変更は、旧版を残して別の試験として記録する。前の凍結した段の試験の変更は §9.1 に登録した2か所だけで、ほかに要ることが分かれば、同じ手順で利用者の了解を得る。

## 付録　凍結の記録（後で記入する）

- 期待する状態の表 `benchmarks/numerics/m4b2/expected_states_m4b2.json` の SHA-256（UTF-8、LF、末尾 LF を含むファイル全体）: `9c97907dab69a00e110414d6862d0d061bc80212000b8f3180f905bcb28108ea`
- 生成のスクリプト `expected_states_m4b2.py` の SHA-256: `58fd1c3d35e2f6e205ae52039cecaf4a0c22ca23a3fe84ea2bde8d533a024d2f`
- 凍結の時の設計 `docs/M4B_DESIGN.md`（第6版）: コミット `54531ce`、SHA-256 `ce853a1084ac7d6021811c7aa1528745f7068d7498c5bf9c4f6b8ca9c67ff191`
- この文書の SHA-256（自己参照の欄を正規化した値）: `304cd057a561b033ec88bbeb1cf2c33cb363e974fb70c4733d995cbf7253921c`。この行の 64 桁だけを ASCII の `0` 64 個に置き換えたファイルのバイト列（UTF-8、LF、末尾 LF あり）に対して計算する。無加工のファイル全体の SHA-256 は、凍結のコミットとその記録に残す。
- R-C の参照 `references_m4b2.json` の SHA-256、生成のスクリプトと SHA-256、環境、所要時間、キーの数（行ごと）、既知の誤答の記録、「値なし」のキー:（参照の担当が、実装の値を見る前に記入する）
- 完了した比較の表は M4b-2 では作らない（`comparison` は常に `not_matched`）。
