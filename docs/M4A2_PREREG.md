# M4a-2 事前登録案　重みの長さと単一損失乗数

2026-11-26、レビュー用、未凍結。基準は `e4e024fb3a4abbca078306c93e39c1f37d1cb845`。`cf69b90` の `src/`、`tests/` はこの版と同じ。設計第3版（`fc2b596`）§9 の U1–U5 を前提とする。実装担当、独立な参照の写像担当、結果を集計する担当が、同じ入力と同じ判定の分母を使うための表である。[設計](M4A2_DESIGN.md)・[M4a の事前登録](M4A_PREREG.md)・[実行行の JSON](../benchmarks/numerics/m4a2/expected_states_m4a2.json)

## 1. 量と入力

表面法線から 0°、固定した初期運動エネルギー \(E_0\)、有限の観測窓 \(W=[a,b]\) を使う。深さ平均には一定長 \(L_w\) を渡し、式(19)の衛星には長さ \(m\) を掛ける。13a の窓面積に追加の減衰 \(c\) を足して作る長さは、この窓と運動量規約に固有である。S2 は同じ窓の平均後の全断面積を使う。M3a と損失積分の数値的な接続を検証する。有限窓から物性の IMFP を決めること、軌道ごとの長さ、多重散乱、斜出射は対象外である。

| 記号 | 定義と単位 |
| --- | --- |
| \(K_j,A_j\) | 深さ平均核（Å⁻¹ eV⁻¹）とその窓積分（Å⁻¹）。列順は `term_13a, term_13b, term_14, term_15_T1, term_15_T2, term_15_T3, term_16_e, term_16_h, no_hole, hole, total`。正孔差分は符号を保つ |
| \(L_w,m,D,e\) | 長さと乗数は Å。段αの分母 \(D=A^\alpha_{13a}+c\)、段βの S2 の分母 \(D=A^\beta_{\rm total}\)。\(c,e\) は Å⁻¹。\(e\) は該当列の推定誤差 |
| \(S,S_{\rm abs},G_D(T)\) | 一次線の積分、その絶対強度の和、損失 \(T\) で領域 \(D\) に移る一次線強度。領域は各ビンと `outside.low, outside.high`。エネルギーは eV、強度単位は \(I\) |
| \(J_0,J_{1,j},J,w\) | 各領域のゼロ損失、衛星3成分、\(J_0+J_{1,\rm total}\) は \(I\)。\(w=mA_{\rm total}\) は無次元の面積比 |
| \(H_j,H^\pm_{j,D}\) | \(maK_j(a)\)（無次元）と \(maK_j(a\pm0)G_D(a\pm0)\)（\(I\)）。報告した \(L_w,m\) を固定した下端の偏感度。`LbK_b` は \(mbK_j(b)\) |
| \(\epsilon,e_f,e_{\rm ret},u_r\) | \(\epsilon=2^{-52}\)。後の3量は各出力と同じ単位。報告した長さ・乗数での推定誤差、乗数区間を伝えた返却誤差、参照の不確かさ |

`RB1–RB8, RB1a–RB4a, CE1–CE6, Z1, N…` は凍結した M4a の状態表の同名行を引く。材料・一次線・出力端は表の名前から binary64 の hex を復元する。行の上書きを既定に重ねる。面積設定は `epsrel=1e-6, area_epsabs=1e-10, panel_limit=400, kernel_epsabs=1e-12, kernel_epsrel=1e-9, kernel_limit=200`。スペクトルは `spectrum_epsabs_fraction=1e-9, max_cells=2000000` を加える。通常 \(E_0=1000\), \(Q_h=1\), `k_total`。CE の実行には元の `area_epsabs=1e-20` を使う。比較の床だけが \(10^{-18}\) である。

`short` は `synthetic_wide`、上記既定、\(W=[5,6], L_w=20\)、スペクトルは `D1000, E_RB1`。CLI は C01/C02 の設定ファイル・一次線形式で同じ有効入力を作る。`derived` は `short` で `weight_length_kind="inverse_13a_window_area", weight_length_A=None`、種は未指定（20 Å）。`external_m` は `short` に `multiplier="external", multiplier_A=7`、`s2` は `short` に `multiplier="inverse_window_area"`。JSON の `c` は `extra_inverse_length_per_A`、`seed_A` は `weight_length_seed_A` の数値。表の `null` の付随パラメータは引数の省略を表し、導出時の必須キーワード `weight_length_A=None` は常に渡す。明示的な無効値を入れる行は §6 に従う。

## 2. 実行行と終了状態の読み方

表の波括弧・添字は各値につき1行に展開する。正本の JSON は各行に `id, entry, inputs, expect, keys, rule, ci` を持つ。順番は `execution_order`、キー集合の展開は §3–§7。`success` は入口が正常に返ること（CLI 成功は終了0・両計算の stdout は `status=complete`）。変異行では正常対照と指定キーの判定がそろって初めて `detected` になる。旧100行の `expect.m4a_state` は対応する旧行の例外型・理由・属性・台帳への外部キーである。記録専用、sdist の skip、未確定も元の数え方に従う。旧 CE5 の凍結後の訂正は M4a §11 のとおり。 `test_driver` と `test_finish` の行では、JSON の `fixture` は入力名であり、注入する数値・順序・観測点の正本は本文 §5–§6 とする。これらの入力は付録の本文ハッシュで固定する。

| ID（JSON の展開単位） | 入口・入力 | 期待状態 | 比較キー、参照、許容 | 実行 |
| --- | --- | --- | --- | --- |
| V1-{入口}-{mode} | 面積・スペクトルの API と両 CLI、short。省略、API の明示 `weight_length`、同値の外部乗数 | 成功 | 基準との対。旧欄全体、新欄の意味。§3.1 のビット・バイト一致 | CI |
| V1-{RB1,RB1a,RB3,RB3a} | 元の API、同値の外部乗数 | 成功 | 同上。Al/Si とδ/ビン一次線 | 全件 |
| V1-old-{旧ID} | M4a の100行の元の入口・入力 | 旧表の状態 | 元の `check_plan`、参照、許容、実効訂正。V01–08 の旧タグ比較もここに含む | 旧 CI 行／全件 |
| V2-{api_area,api_spectrum}-{external,inverse_window_area,weight_length} | short の窓13a重み、\(c=1/32\)、外部乗数は 7 Å | 対が成功 | 同じ規約で段βを再実行、§3.2 | CI |
| V2-real、V2-real-weight_length | 窓13a重み、\(c=0\)。順に RB1 のスペクトル＋S2、RB1a の面積＋乗数 `weight_length` | 対が成功 | §3.2。RB1a は固定時の誤差を対照と比べ、導出側の \(H_j,w\) の返却誤差を §4.4 でビット照合 | 全件 |
| V2-area-only | 面積、short、\(m=7\) と省略 | 両方成功 | 11面積、誤差、予算、各内訳をビット比較 | CI |
| V3-{RB行}、V4-{RB行} | 元の12行の面積／スペクトル API。\(L_w=L_{\rm ref}\)。V3 は外部 \(m=2L_{\rm ref}\)、V4 は S2 | 成功 | 該当する R-B の全キー、ビン密度、§4。V4 では固定時と区間込みを別々に判定 | 各 RB1 のみ CI |
| V5-{RB/CE行}-{z,p}-{10,20,40} | 18行それぞれの段α入口。z は \(c=0\)、p は \(1/32\)、末尾は種（Å） | 種の段が成功 | `term_13a` の参照、\(L_w\) と区間、§4.5。段βの成功をこの行の条件にしない | RB1a-z-20 のみ CI |
| V5-{RB/CE行}-{z,p}-intersection | 同一の参照行・\(c\) の上記3種 | 3実行が成功 | 閉区間の共通部分が非空、各種個別にも照合 | 全件 |
| V3-layer、V3-tails | `integrate_window`、定数核／強度差の大きいδ線 | 成功 | 全領域・密度・感度／外側2領域の直接積分、§5 | CI |
| V4-stages、V4-error | 試験入口。2段の実駆動／非対称区間で仕上げ | 成功 | 各段の \(m\)、残したパネルの値と5項／補正の式、§4.4、§5 | CI |
| V6-wire、V6-{inside,equal,outside} | `_run` の核の実引数／α・βの境への注入 | wire、内側、等号は成功。外側は収束例外 `weight_length_inconsistent`、production 1（βの最終合格段） | 長さ、βの列、整合差と許容を有理数で比較、§5 | CI |
| V7-{a1,a2,ext} | 層入口。定数核、窓下端 1 または 2。ext は \(m=4\) | 成功 | \(w,H,H^\pm,\mathrm{LbK}_b\)、収支を解析値と比較、§5 | CI |
| V8-{Al,Si}-{default,historical} | 面積、窓の13a、\(c=0\)、種20または当時の種、外部 \(m\) | 成功 | 保存した \(L_w,w\)、§5 の固定幅 | 全件 |
| V9-* | 面積・スペクトル・試験入口。型、順序、構造的ゼロ、段ごとの分母、細分・上流の失敗 | §6 の各行 | 例外型・理由・phase・stage、実呼出し・分割、厳密な逆数と区間 | CI |
| V10-* | 7規約×両 CLI と対応する API、再生、型、照合と誤った旗 | §7 の各行 | 記録・単位・警告・状態、同一実行の ID、再生のビット一致 | CI |
| V11-M01–M26 | 一変異ずつ §8 の入口で正常対照と対に実行 | 指定の失敗を検出 | §8 の行・キーに対する同じ規則 | 全件 |

## 3. 互換と同条件再実行

### 3.1 V1 の対

同じインタプリタ、OS、NumPy/SciPy、浮動小数設定、環境変数、入力バイト列を使う。`git rev-parse e4e024f^{commit}` の完全 ID を 文頭と照合し、`git archive --format=zip <完全ID> src/xpsloss` を一時展開する。基準と試験対象を別の子プロセスで1回ずつ起動し、それぞれの `src` を先頭に置き、`xpsloss.__file__` の絶対パスがその配下であることを確認する。各子の中で複数行を実行してよい。基準側には旧来の引数だけを渡し、変更版で各 mode を与える。CLI に同長規約を明示する旗はない。同値の外部指定は API で `multiplier="external", multiplier_A=L_w`、CLI で `--multiplier-A L_w`。Git checkout で基準が取得できなければ失敗、展開した sdist だけはこの対を理由つきで skip とする。

- 面積の旧欄は `electron_energy_eV, loss_window_eV, weight_length_A, epsilon_momentum, hole_charge, theta_deg, window_area, window_area_error, window_area_budget, single_loss_weight, single_loss_weight_error, inverse_window_area_A, inverse_window_area_interval_A, inverse_reason, cutoff_sensitivity, cutoff_sensitivity_error, diagnostics, warnings`。
- スペクトルの追加の旧欄は `output_edges_eV, primary, primary_integral, primary_abs_integral`、5量 `zero_loss, single_loss_total, single_loss_no_hole, single_loss_hole, spectrum` のそれぞれに `（無接尾）, _error, _per_eV, _per_eV_error`、`cutoff_sensitivity_bins, cutoff_sensitivity_bins_error, outside`。一次線も全 dataclass 欄、両側の外も全子キーを比べる。
- `diagnostics` は `resolution_basis, resolution_basis_reason, candidates, remote_root_checks, required_panels, final_panels, stages, gradings, tolerance_bisected_panels, verification, error_budget, balance, kernel_calls, effective_settings, comparison, warnings, LbK_b`、スペクトルの `primary_support, initial_cells, cells`。設計 §2.4 の1–2行は、API と CLI 記録の `diagnostics` について全階層のキー集合が基準と完全に一致する。`multiplier`、`denominator`、`weight_length` の診断ブロックは追加しない。`warnings` は上流警告の出現回数・点・順を含む。辞書は旧キーごとの型と値、列は型・順・長さ、配列は dtype・形・要素ビット、float は binary64 のビット列（符号つきゼロを含む）、整数と真偽は型も区別する。
- CLI は stdout の `status, model_id, output_dir, files, warning_count`、stderr と終了、作成ファイル集合を比べる。`window_area.csv, single_loss_spectrum.csv` はヘッダ・改行を含めバイト一致。`run.json` は `schema_version, status, command, created_utc, elapsed_seconds, input_sha256, input_file, package_version, model_id, dependent_model_ids, runtime, config, units, constants, files, warnings, spectrum_input, spectrum_definition, effective_settings, effective_settings_sources, primary, output_edges_eV, results, diagnostics, comparison` の型と全子孫。時刻は UTC 表現を parse、経過秒は有限かつ非負、入力・出力・一次線ファイルの一時パスは各実行の実在する対応物に展開して basename・役割・入力ハッシュを照合する。これらだけを値比較から除く。
- 追加を認める位置は設計 §3–4 で定めた新しい結果欄、記録・単位等（1–2行の `diagnostics` は基準と完全一致）。旧警告列・`effective_settings` のキー集合と来歴・旧照合は完全一致。追加キーを消してから新規欄の試験まで通したことにはしない。§7 で個別に照合する。凍結した旧100行は行別に既存の検査を呼び、再集計時も失敗・未確定を残す。上流をそのまま送出した場合の `vars()` 全体比較（`tests/m4a_states.py` の `upstream_equal`）は上流例外に属性を付加すると失敗する。N29 の `check_n29` は包んだ例外の内部属性を列挙して比べるため、外側の `phase` は null とする。V1-old-N25 は元の検査に加え、変更版で乗数を省略した場合と同値を外部指定した場合の両方で同じ面積計算を失敗させる。どちらも `loss_budget_not_met` の旧属性・診断を保ち、例外に `phase` と `stage` が存在し、両方とも null であることを読む。

### 3.2 V2 の段β

導出側の確定した `weight_length_A` を対照の `external_constant` に渡す。種と \(c\) は外し、窓・一次線・端・材料・全実効設定・電荷・運動量をそろえる。**乗数の規約は同じ**にする。S2 は対照も S2 を計算する。

1. §3.1 の旧結果のうち、入力の数値と一次線、面積11列・誤差・予算、逆面積・区間・理由、\(w,H\)、5量、密度、両側感度、左右の外と各誤差を、各入口が返す範囲でビット比較する。返却誤差のうち手順2の同長対には同手順の規則を使う。段βの `LbK_b, resolution_basis, resolution_basis_reason, candidates, remote_root_checks, required_panels, final_panels, stages, gradings, tolerance_bisected_panels, verification.{rho_max,epsilon_max}, initial_cells, cells, primary_support, effective_settings, balance.*.{difference,rounding_bound}`、`kernel_calls.{calls,points_sent_per_stage,points_sent_total,distinct_points,batch_sizes}`、`error_budget.outputs, error_budget.max_ratio` も同じ。S2 は段βの分母の値・誤差・閾値・確定段も同じ。
2. `weight_length` の対だけは、導出側に乗数区間があり、対照側の区間は null になる。両側で**固定時**の誤差 \(e_f\) を直接読み、ビット比較する。導出側の区間補正対象の返却誤差（密度・合計・\(w\) を含む）は §4.4 で検算する。面積と \(J_0\) の返却誤差はこの対でもビット一致。外部乗数同士、S2 同士では返却誤差もビット一致。対照が設計 §2.4 の1行目になる同長対では `diagnostics["multiplier"]` がないため、導出側だけで区間と `rho` を照合する。`model_id`、規約・区間・来歴、段α、警告、照合は各実行について §7 を使う。
3. V2-area-only は `window_area, window_area_error, window_area_budget` の11列、`error_budget.outputs["window_area."+j]` の `interpolation, verification, kernel_point_error, edge, rounding, estimate, budget, ratio, flag` をビット比較する。\(w,H,\mathrm{LbK}_b\) は新しい \(m\) で計算した値を §5 の解析例でも確かめる。

## 4. 独立な参照への写像と判定

[凍結した R-B](../benchmarks/numerics/m4a/references_rb_m4a.json) の `rows[id].keys[key]` を、同じ行の `conditions` で読む。R-B は M3a の点計算を共有し、損失積分と一次線・領域の積分を別経路で求めている。CE はその経路の面積列である。V2、V8、行内の恒等式は独立な数値比較の件数に入れない。

### 4.1 独立担当が作る表

本文と状態表の凍結後、実装・試験の観測値を見ない別担当が作る。この下書きには写像済みの数値を収めない。1レコードの主キーは `(row_id, check_kind, key)`。`check_kind` は `external, fixed_m, ratio, seed_area, inverse_interval` 等で、一つの出力の二つの比較を上書きしない。各レコードに次を保存する。

| フィールド | 内容 |
| --- | --- |
| `source_sha256, source_row, source_keys, transform, unit` | 元ファイルのハッシュ、分子・分母・一次線のキー、演算の種類、単位。領域・±・列の添字を含む |
| `value, u_r, d_abs, d_rel, interval` | 中心、非負の半幅、絶対・相対の床、必要な場合の外向き両端。\(r=0\) の `d_rel` は null。逆数の区間に端の架空の誤差を付けない |
| `status, reason, parameter` | `value / parametric / no_reference / reference_path_failure` と理由。条件付き比較の `parameter` には \(\hat m\) の記号を残す |
| メタデータ | スキーマ版、入力と生成器のハッシュ、精度と丸めの規則、行・種類別の予定キー数、生成時の件数。比較前にファイルの SHA-256 を別記する |

元の十進文字列と binary64 入力はそれぞれ厳密な有理数として読む。積・商・和・ビン幅を有理数で処理し、保存は十進の有効数字 **40桁**（0 は `"0"`）。中心は最近接・同率偶数、端は外向き、\(u_r\) は元の区間を覆う半幅に中心の十進丸めを含めて切り上げる。丸めが厳密な項に架空の幅を足さない。比較の床は下方丸めまたは厳密な有理式も記録してその式で再現する。\(d_{\rm rel}\) の上限判定は丸める前と保存後の両方で確認する。キーごとに必要な参照の入手可否を判定する。

V4 の `fixed_m` だけは凍結時に \(\hat m\) が未定である。基準長、\(q\) に掛かる分子の中心・幅、加える \(J_0\) の中心・幅、床の係数と演算順を凍結する（定数項と比例項の `value,u_r`、`d_abs,d_rel` の評価式）。比較時に**一度だけ**報告した binary64 の \(\hat m\) を代入し、その hex と評価後の数値4欄を実行台帳に残す。独立担当の表の書き換えは生じない。`ratio` は全面的に参照だけで事前に数値化できる。

### 4.2 V3 と V4①　固定した乗数での比較

\(L=L_{\rm ref}\) は同じ行の `conditions.weight_length_A`。V3 は指定した外部 \(m\)、V4① は固定した \(\hat m\) を用い、\(q=m/L>0\)。元キーの中心を \(r_i\)、閉区間を \(I_i=[r_i-u_i,r_i+u_i]\) とする。まず下表で中心と区間を同じ式で写す。左右の外は元の同じ側の直接積分から取る。和の幅は \(u_0+qu_1\) と保存時の丸めである。

| キー | 中心の式 | \(d_{\rm abs}\) |
| --- | --- | --- |
| `window_area.j` | \(r_j\) | \(B_A^0=\) 元の `area_epsabs`（CE は \(10^{-18}\)） |
| `single_loss_weight` | \(qr_w\) | \(mB_A^0\) |
| `cutoff_sensitivity.j` | \(qr_H\) | \(ma\,10^{-12}\) |
| `D.zero_loss` | \(r_0\) | \(\eta S_{\rm abs}\) |
| `D.single_loss_{total,no_hole,hole}` | \(qr_1\) | \(\eta S_{\rm abs}\) |
| `D.spectrum` | \(r_0+qr_{1,\rm total}\) | \(\eta S_{\rm abs}\) |
| `D.cutoff_sensitivity.{plus,minus}.{total,no_hole,hole}` | \(qr_\pm\) | \(ma\,10^{-12}S_{\rm abs}\) |
| ビンの5密度 `bin[i].量_per_eV` | 上表の写像済みの量を同じビンの厳密な幅 \(\Delta E>0\) で割る。区間も割る | 対応する床を \(\Delta E\) で割る |
| `inverse_window_area_interval_A.{lower,upper}` | 参照の `window_area.total` の区間を反転。元の面積・誤差で丸めも検算 | §4.5 の区間判定 |

\(\eta=\) 実行の `spectrum_epsabs_fraction`。個々の比例キーでは \(u'=qu\) なので、十進丸めを除けば相対幅は元のまま。和・比・相殺を含め、適格性と \(d_{\rm rel}\) は必ず写像後に再計算する。V3-RB1 の `bin[6].spectrum` と V3-RB3 の `outside.high.spectrum` は、\(q=2\) のときもそれぞれ 1 と約 0.35（衛星0）。両者の `zero_loss` と `single_loss_total` も個別に比較する。

### 4.3 V4②　S2 の比の比較

段βで使った参照条件の \(A=r_{A,\rm total}\)、区間 \([A_-,A_+]\) を **同じ凍結行** から取り、\(A_->0\) を先に検査する。\(m_{\rm ref}=1/A\)、中心の \(q=1/(LA)\)。床の \(m\) は \(m_{\rm ref}\)。\(c\) はこの分母に加えない。外部 \(L=L_w\) の行だけをこの写像に割り付ける。未知の \(L_w\) で凍結 R-B の全核を代用することはできない。

\(m\) に比例するキーの分子を \(N=r_1/L\)、その区間を \([N_-,N_+]=I_1/L\) とする。比の中心は \(R=N/A\)。区間の下端と上端は次の4数の最小と最大である。

\[
 \{N_-/A_-,\ N_-/A_+,\ N_+/A_-,\ N_+/A_+\}.
\]

符号つき `hole`、0 をまたぐ感度、片側の外も同じ演算を使う。非負を仮定した端の選択、独立誤差を仮定した二乗和は用いない。`spectrum` はこの衛星全量の区間に \(I_0\) を端ごとに加える。面積・\(J_0\) は不変、密度は和の後で幅で割る。中心から両端までの大きい距離を \(u_r\) にする。参照分子と分母は別の積分値でも、同じ上流誤差を共有し得る。S2 の \(w\) だけは同じ面積が消去されるため、②では標的 1 との丸めの恒等式として計上する。①の固定時の \(w\) の比較には \(qr_w\) と加算前の \(e_f\) を使う。

### 4.4 判定の優先順と返却誤差

1. 必要な元キーが欠ける、元キーが `no_reference`、分母の参照区間が 0 を含むか負なら、依存するキーだけ `no_reference` とする。元の参照経路自体の失敗は `reference_path_failure`。いずれも予定の行・キーを削除しない。\(r\ne0\) なら \(d_{\rm rel}=\max(8u_r/|r|,\varepsilon_{\rm rel})\)。\(10^{-4}\) 超は `no_reference`（上限に切り詰めない）。\(r=0\) は \(u_r\le d_{\rm abs}/4\) が必要で、正の床で \(u_r>d_{\rm abs}/4\) なら参照不足。床0で \(u_r>0\) も同じ扱いとし、確定した厳密な0に対して \(y\ne0\) なら不合格。
2. 有効な参照で \(\delta=|y-r|\)、\(d=d_{\rm abs}+d_{\rm rel}|r|\)（\(r=0\) は絶対床のみ）。**精度**は \(\delta+u_r\le d\)。**推定の被覆**は \(E=e_y+4\epsilon\max(|y|,|r|)\) として \(\delta+u_r\le E\) は被覆、\(\delta-u_r>E\) は過小、それ以外は部分的な重なり。M4a §7.2–7.3 と同じく精度不合格、非有限値・負の誤差、過小推定、必要な予算の違反は **不合格**。精度と必要な予算を満たし、被覆なら **合格**、部分的な重なりなら **未確定**。精度と被覆それぞれの判定も保存する。
3. 予算対象は全11面積と各領域の \(J_0,J_{1,\rm total},J_{1,\rm no\_hole},J_{1,\rm hole}\)。面積は \(B=\max(\text{area_epsabs},\varepsilon_{\rm rel}|y|)\)、領域は \(B=\max(\eta S_{\rm abs},\varepsilon_{\rm rel}|y|)\)。V4① は試験入口から直接取り出した \(e_f\) で §4.2 の精度・被覆と \(e_f\le B\) を判定する。派生キーに個別予算を作らない。V4② は公開の \(e_{\rm ret}\) で §4.3 の精度・被覆を判定する。返却誤差に固定時の予算を課さない。
4. V4③・V2 では、乗数区間を持つ側の返却誤差を次の順でビット照合する。各式の引数は記録された binary64 を厳密な有理数に戻す。`up` は厳密値以上の最小の binary64、`down` は厳密値以下の最大の binary64、`RN` は最近接・同率偶数。非負の厳密なゼロの返却誤差は +0 とする。

   | 比べる返却キー | 期待する binary64 |
   | --- | --- |
   | `diagnostics["multiplier"]["rho"]`（乗数区間あり） | \(\rho=\mathrm{up}(\max(\hat m-m_-,m_+-\hat m)/\hat m)\) |
   | 各ビンの `single_loss_{total,no_hole,hole}_error`、`outside.{low,high}.single_loss_{total,no_hole,hole}_error`、11列の `cutoff_sensitivity_error[j]`、`cutoff_sensitivity_bins_error[sign][component]` と `outside.{low,high}.cutoff_sensitivity_error[sign][component]` | 各 \(\hat y,e_f\) と上で丸めた \(\rho\) を使い、\(e_{\rm ret}=\mathrm{up}(e_f+\rho(\lvert\hat y\rvert+e_f))\)。積と和は厳密に求め、最後に1回だけ上向きに丸める。`sign` は plus/minus、`component` は total/no_hole/hole |
   | `single_loss_weight_error`（窓13aの重み＋乗数 `weight_length`） | \(\hat y=w\)、同じ \(m\) での固定時の誤差 \(e_f\)（手順5）を上と同じ式で補正 |
   | 各ビンの `spectrum_error`、`outside.{low,high}.spectrum_error`（乗数区間あり） | \(\mathrm{up}(e_{1,\rm ret}+e_0+\epsilon\lvert J\rvert)\)。報告する衛星全量の補正済み誤差・ゼロ損失誤差・報告した合計強度を使い、厳密な和を1回だけ上向きに丸める |
   | ビンの `{zero_loss,single_loss_total,single_loss_no_hole,single_loss_hole,spectrum}_per_eV_error` | 対応するビンの最終誤差と値を使い、\(\mathrm{RN}((e_{\rm bin}+\epsilon\lvert y_{\rm bin}\rvert)/\Delta E)\)。\(\Delta E\) は両端の厳密な差。乗数補正は各積分量で1回だけ行う |
   | `window_area_error[j]`、各ビンの `zero_loss_error` と `outside.{low,high}.zero_loss_error` | 固定時の同じキーのビット列 |

   乗数に比例する各キーでは、厳密な \([\hat y-e_f,\hat y+e_f]\,[m_-/\hat m,m_+/\hat m]\) の四隅がすべて \([\hat y-e_{\rm ret},\hat y+e_{\rm ret}]\) に入ることも検査する。丸めの照合と包含のいずれにも ulp の許容幅を足さない。
5. S2 の `single_loss_weight_error` は乗数補正と別に扱う。\(\hat m=\mathrm{RN}(1/A^\beta_{\rm total})\)、報告 \(w=\mathrm{RN}(\hat m A^\beta_{\rm total})\)、誤差 \(=\mathrm{up}(|w-1|)\) をビットで検算する。衛星全量を全ビンと左右の外で足した値と \(S\) の差は、`balance.single_loss_total.rounding_bound` と \(S|\hat m A^\beta_{\rm total}-1|\) の和以下。外部から与えた乗数または外部重みとの同長規約での \(w\) の誤差は M4a と同じ演算順の \(me_A+\epsilon|w|\)、窓13aの長さを乗数に使うときはこの \(e_f\) に手順4の補正をする。固定 \(L_w\) の核の推定だけが含まれる。
6. 単位は `(row_id,check_kind,key)`。精度、被覆、予算、丸め、恒等式、変異検出を別列に持つ。`pass, fail, unresolved, no_reference, not_run`、予期しない実装未収束、上流失敗、参照経路失敗を分け、母数・件数を行と種類別に示す。区間の部分重なりは `unresolved` の内訳であり参照不足ではない。実装が止まって読めなかったキーは `not_run`、行には実際の例外を残す。V3・V4 の `_per_eV` は同じ参照のビン積分を幅で割るので、合否を保存し、独立な数値比較の件数から分けて集計する。各写像の RB 12行の予定1,308キーの内訳は、直接参照に由来する1,024、逆数区間の端24、派生密度260（52ビン×5量）。V4①・②はこの内訳を別々に報告し、②の \(w\) は手順5の恒等式に計上する。V5 の3種比較で一つが未確定なら、それを合格へ繰り上げない。

### 4.5 V5 と逆数区間

1. R-B 系12行、CE 系6行それぞれについて \(A^\alpha_{13a}\) を `rows[id].keys["window_area.term_13a"]` の \(r,u_r\) と **面積の単位で** §4.4 により比較する（\(e_y=e^\alpha\)）。\(c\) はこの面積比較には入らない。
2. \(\mathrm{Fr}\) を binary64 からの厳密変換とする。\(D=\mathrm{Fr}(A^\alpha)+\mathrm{Fr}(c)\)、\(e=\mathrm{Fr}(e^\alpha)\)。\(D-e>0,\ e\le\mathrm{Fr}(\varepsilon_{\rm rel})|D|\) を比べ、長さを \(\mathrm{RN}(1/D)\)、区間を \([\mathrm{down}(1/(D+e)),\mathrm{up}(1/(D-e))]\) とビット照合する。最近接は同率偶数。**和を binary64 で先に丸めない。** \(c\) を 0 にして全面積で同じ検算を行うと、V3/V4 の `inverse_window_area_A` とその区間になる。S2 の `multiplier_interval_A` はその区間と同じ。
3. 参照の下端 \(r-u_r+c>0\) のとき、参照区間 \(I=[1/(r+u_r+c),1/(r-u_r+c)]\) を厳密なまま報告区間 \(P\) と比べる。\(I\subseteq P\) なら合格、\(P\cap I=\varnothing\) なら不合格。それ以外（端だけが接する場合を含む）は未確定。手順2が不一致なら包含しても不合格。面積の許容を反転前の不確かさに足さない。種10、20、40 Å は各々で手順1–3を判定し、同一条件で \(\max L_-\le\min L_+\) も求める。
4. 登録は18参照行×2通りの \(c\)×（3種＋交差1行）＝144行。13a 面積の独立な条件は材料（ギャップと振動子）×窓×運動量規約×\(E_0\) で集計し、全件では **11条件**。RB1・RB1a・RB5・RB6・RB7 は同じ Al の 13a 面積になる。RB2/2a、RB3/3a、RB4/4a も各1条件で、\(E_0\) が異なる RB8 を加えて RB は5条件、CE は6条件。電荷・長さ・一次線・\(c\)・種の変更による照合件数は面積の独立条件数に加えない。集計時も `reference_inventory.v5_area_groups` の同一条件内にある失敗・未確定を個別に残す。

## 5. 配線、感度、再現に使う入力

合成核は binary64 の二進分数で11列を作る。特記しない項と点の誤差は0、`no_hole`、`hole`、`total` は M4a の列和と \(Q_h=1\) で作る。注入する入口は内部の観測用でよい。実装担当が固定するのは実引数、観測量、終状態であり、新しい公開 API 名は必要ない。

| 行 | 入力・手順 | 判定する量と幅 |
| --- | --- | --- |
| V3-layer、V7-a1/a2/ext | `constant` は \(K_{13a}=1/2,\ K_{14}=-1/4\)、したがって \(K_{\rm total}=1/4\)。\(L_w=8\)、\(b=4\)、既定 \(a=2\)、\(p=10,s=1\) のδ線、端 [6,8,9,11]。V3-layer と V7-ext は \(m=4\)。層への \(m\) は既定つきキーワードまたは実際の駆動経路の試験入口で指定 | \(A_j=(b-a)K_j,\ w=mA_{\rm total}, H_j=maK_j,\ \mathrm{LbK}_{b,j}=mbK_j\)。各領域の衛星は \(mK_j\) × \(W\cap(p-D)\) の長さ、\(J_0\) は \(p\) の所属領域、\(J=J_0+J_{1,\rm total}\)。各領域を直接比べ、密度はビン幅で割る。解析値を \(r\) として \(|y-r|\le8\epsilon |r|\)、0 は厳密 |
| V7-a1/a2 | S2 をそれぞれの窓で解き直す。\(T=a+0\) なら \(p-a\) のすぐ下側、\(T=a-0\) ならすぐ上側のビンに入る。\(a=1\) で境は9、\(a=2\) で境は8 | 両窓で \(w\simeq1\)、\(H_{\rm total}=a/(4-a)\)。\(H^\pm_{j,D}=maK_j\) が該当する片側のビンだけに入り、反対側は0。収支の `difference, rounding_bound` は M4a の `_balance` に \(m\) を代入した式で照合（§5 下段） |
| V3-tails | 定数核、\(W=[2,4],m=4,L_w=8\)、δ位置 [10,100,200]、強度 \([1,2^{60},2]\)、ビン [90,110] | 両外の \(J_0=(1,2)\)、\(J_{1,\rm total}=(2,4)\) を直接検査。幅は各側 \(8\epsilon|r|\)。全領域から一方を逆算すると、binary64 の全強度から小さい尾が失われる |
| V4-stages | 実際の段βを \(W=[2,4]\)、初期パネル [2,3],[3,4]、\(p=10,s=1\)、端 [6,7,8,9,11]、\(K_{13a}=1\) で駆動する。試験注入で初段の面積集計の `total` だけを3（誤差0）にし、非分母列の局所誤差で左1枚だけ二分させる。次段は注入を外し真の面積2。両段とも分母は確定 | 各段の `fold` 実引数は \(\mathrm{RN}(1/3),\mathrm{RN}(1/A_{\rm final})\) で異なる。未分割の [3,4] が寄与する [6,7) の衛星は最終 \(m\) × 1。2段目に得た全パネルを同じ最終 \(m\) で個別に畳む固定パネル対照と、領域の値・誤差5項・局所の配分・予算をビット比較。分割回数、存続パネル、呼出時の \(m\) を記録 |
| V4-error | `test_finish` に \(D=2,e=1,\varepsilon_{\rm rel}=1/2,\hat m=1/2\)、外向き区間を与える。符号つき衛星 \(\hat y=-2,e_f=1/8\)、別キー \(\hat y=0,e_f=1/8\)、\(J_0=3,e_0=1/16\)、幅2のビン、同じ値の外側・片側感度を用意する。仕上げのみを試す | \(\rho=1\)、補正後はそれぞれ \(9/4,1/4\)、\(e_0\) は元のまま。`spectrum` の衛星全量にも非零の値と誤差を与え、合計・密度・\(w\) を §4.4 で検算する。全量2、誤差1/8なら \(J=5\)。さらに同じ区間で衛星の別のビンに \(\hat y=1,e_f=2^{-54}\) を置く。厳密な誤差 \(1+2^{-53}\) の返却ビット列は `0x1.0000000000001p+0`（上向き）。丸めの向きをこのビンで判別する |
| V6-wire | `_run` の `_call_kernel` を長さ依存の合成核で置き換える。\(W=[2,4],L_{\rm seed}=20,c=0,m=7,\ K_{13a}=1/8,\ K_{13b}=L_{\rm arg}/64,\ K_{14}=-1/32\)。能動材料、窓・運動量・点・設定と段ごとの `weight_length_A` を spy で保存する | αは20、βは \(\mathrm{RN}(1/\hat A_{13a}^\alpha)\)（解析値4 Å）。13bの返却面積は \(2L_w/64\)（種なら \(40/64\)）。解析値は \(8\epsilon|r|\)、実引数と回数は厳密。上位 `kernel_calls` はβのみ。αの13a以外の列だけに予算超過の点誤差を入れ、αがそれで細分しないことも同時に見る。βでは注入を外す |
| V6-inside/equal/outside | 試験入口で最終の \(A^\alpha=1/4,\ e^\alpha=e^\beta=2^{-24}\)。\(x=1/4+2^{-23}\) とし、βの面積を順に `nextafter(x,0), x, nextafter(x,+inf)` とする。窓 [2,4] の導出重み、\(c=0\)、乗数 `weight_length`。αは初段で成功。βは初期1枚・上限2枚とし、outside だけ初段の `term_13b` の面積0に検証誤差 \(2\,\text{area_epsabs}\) を与えて1回二分する。次段でその誤差を0にする。各行で最終段の全列の予算は通る | 入力の binary64 を有理数に戻し \(\lvert A^\alpha-A^\beta\rvert\le e^\alpha+e^\beta\) の真偽を比較する。外側では §2 の例外の `stage` が最後に予算を通ったβの整数段（この入力では1）と一致し、ほかは成功。診断の両面積・誤差・差・許容も検算 |
| V8-Al、V8-Si（各2種） | ELF DB v1.0 の `benchmarks/inputs/{al,si}_elfdb_v1.0.bulk.json`、1000 eV、法線、電荷1、`k_total`、\(c=0\)、面積設定6個は §1 既定。Al は [1,50], \(m=19.4\)、Si は [3,50], \(m=24.6\)。当時の種は各 \(m\)、既定種は20 Å | Al の \((L_w,w)=(27.0571259020694,1.0861615888385654)\)、Si は \((25.549842388377478,1.6938314906967444)\)。各量で \(|y-r|\le4\times10^{-6}|r|\)。`result_C-{Al,Si}.json` の `B_length` と `checks.B.M4a` × \(m\)。固定深さを数値平均した `fine.B` の値はここに使わない |
| V10-coincidence | `test_finish` に外部 \(L_w=20\)、S2、\(A_{\rm total}=\mathrm{RN}(1/20)\)、\(e=0\) を渡す。確定する乗数は20 | §7 の行4として模型 ID・警告・照合・乗数区間を照合する |

収支の直接検査は領域順の値を \(x_D\)、M4a での各領域の `edge+rounding` の和を \(R\) とする。衛星では \(t=\mathrm{RN}(Sm A_j)\)、\(R\) に \(Sm\) × 当該面積の丸めを加える。合計では \(t=S(1+w)\)、面積は total。ゼロ損失では \(t=S\)。`difference` は `fsum` で \(\sum x_D-t\)、`rounding_bound` は \(R+8\epsilon(\sum|x_D|+|t|)\)。基準版と同じ演算順で最終の報告 \(m\) を用いて照合する。収支が通っても左右の外を合格にする根拠には使わない。

V8 の幅は当時の面積精度 \(10^{-6}\) を単位に、種によるαの細分差とβの積分差に各2単位を配した再現用の許容である。保存したαの相対推定は Al \(5.81\times10^{-8}\)、Si \(8.33\times10^{-10}\)、βの全量は順に \(5.15\times10^{-8}\)、\(4.50\times10^{-10}\)。種20 Å と当時の種で細分が変わっても同じ幅で照合する。RB3 の Si は窓下端5 eV なので V8-Si の独立参照に流用しない。

## 6. V9 の拒否、保留、上限

以下の I は `SingleLossInputError`、C は `SingleLossConvergenceError`。この2種類の例外について、設計 §2.4 の組み合わせと失敗時点から `phase` を決める。前検査・段αを飛ばした構造的ゼロでは null。1–2行の既定経路は計算中の例外でも null。3–7行の計算中は、種の計算なら `"weight_length_seed"`（α）、本計算なら `"production"`（β）とする。

`stage` は各 phase 内で0始まり。M4a から継承する理由は、どの経路でも元の例外属性を保ち、`phase` だけを加える。`kernel_not_converged` は核を呼んだ整数段、`loss_budget_not_met`、`kernel_error_exceeds_budget`、`panel_width_unresolved`、`cell_limit_reached` の `stage` は null。新しい理由は、計算前なら null、反復中ならその整数段、`weight_length_inconsistent` は条件を満たした段βの最後の整数段とする。上流から直接返る例外・素の `ValueError` の属性は基準版どおり。`kernel_evaluations=0` は実際の核呼出しと送信点数がともに0であること。新しい理由のメッセージには対応する引数名を含める。

| V9 の ID 接尾辞 | 入口と入力（1値1行） | 期待する終了状態とキー |
| --- | --- | --- |
| {weight_length_kind,multiplier}-{unknown,np_bool,array2,none} | short、該当引数に列挙外文字列、`np.bool_(True)`、`np.array([0,1])`、None | 前検査 I、`invalid_weight_length_kind` または `invalid_multiplier`。文字列比較前の型検査 |
| {extra_inverse_length_per_A,weight_length_seed_A,multiplier_A}-{token} | c と種は derived、乗数は external_m。token は `"2", True, np.bool_(True), NaN, +inf, -1, 2**2000`（binary64 変換が overflow）、種と乗数には 0 も | 前検査 I、順に `invalid_extra_inverse_length, invalid_weight_length_seed, invalid_multiplier_length` |
| old-length-{token} | 外部重みの長さに None、`"2"`、両種の真偽、NaN、+inf、0、−1 | 前検査、素の `ValueError`。基準コードと型・文字列一致 |
| derived-length / external-c / external-seed | derived に長さ20／外部重みに \(c=1\)／外部重みに種10 | 前検査 I、`invalid_weight_length / extra_inverse_length_without_derived_weight / weight_length_seed_without_derived_weight` |
| m-missing / m-extra / s2-extra | `external` に乗数 None／既定に乗数20／S2 に乗数20 | 前検査 I `invalid_multiplier_length` |
| order-kind / order-old / order-window | kind 不正と長さ負／外部長負と乗数規約不正／導出指定、窓反転と乗数規約不正 | 前検査、順に I `invalid_weight_length_kind`、旧版の素の ValueError、I `invalid_window` |
| api-zero-c / pre-N30 | 外部重みに明示 \(c=0\)／公開スペクトル API に旧 N30 の入力で導出を指定 | 成功して旧欄一致／前検査 I `cell_limit_exceeded`（初期89セル、上限88） |
| pre-N35a | 公開面積 API。旧 N24 の材料 N1、\(W=[2.1,10]\) eV、`panel_limit=400`、重みを窓13aにする | 前検査 I `resolution_exceeds_panel_limit`、`required_panels=401` |
| pre-N32a | 公開面積 API。derived の窓を \([1,\mathrm{nextafter}(2,+\infty)]\) eV にする。窓全体の節点は分離でき、初期分割が 2 eV に境を作る | 前検査 I `panel_width_unresolved`。例外の `panel_eV` は \([2,\mathrm{nextafter}(2,+\infty)]\)、`panel_width_ulp=1` |
| zero-w / zero-both / zero-s2 / zero-c-s2 | 構造的ゼロ Z1。順に導出 \(c=0\) 同長、導出 \(c=0\) S2、外部重み S2、導出 \(c=3\) S2 | 前検査 I。順に `weight_length_undefined, weight_length_undefined, multiplier_undefined, multiplier_undefined`。最後は長さが確定可能でも全量は構造的に0 |
| zero-c-{weight_length,external} / zero-overflow | Z1、導出 \(c=3\) 同長または外部7／\(c=2^{-1074}\) 同長 | 成功、\(L_w=\mathrm{RN}(1/3)\)、非同一の外向き端、種を飛ばす `structural_zero`、αの呼出・面積・誤差0、整合 `not_applicable`、\(J=J_0,w=0\)／C `weight_length_not_representable`、phase=null、stage=null、核呼出し0（段αを飛ばした構造的ゼロ） |
| zero-primary | derived、\(c=1/32\)、S2、一次線 `D1000zero`、活動材料 | 成功。面積と分母は有効、\(w\simeq1\)、5量は全領域で0。ゼロ予算・ゼロ誤差の比は null。観測の比は 0/0 |
| signed | `_run` で能動材料に §5 の定数核を注入し、同じ窓 [2,4]、同じ一次線・出力端、導出重み、\(c=1/8\)、S2 を使う | 成功。\(L_w\simeq8/9,\ m\simeq2,\ w\simeq1,\ A_{\rm hole}=-1/2\)。衛星の負のビンを領域ごとの解析積分と比較する。解析値との幅は §5 の定数核と同じ |
| {a,b}-negative | 下記の段入口、初段で \(D=-1,e=2^{-24}\)、実際の核呼出し後に集計値へ注入 | I `weight_length_undefined / multiplier_undefined`、α/β、stage=0。\(D+e<0\) を確認 |
| {a,b}-positive-limit / {a,b}-zero-limit | 正の未確定分母で初段から分割不能／活動材料で \(D=e=0\) が二段続き、次の分割が不能 | C `weight_length_unresolved / multiplier_unresolved`、α/β、stage=0／1。既存の最悪キー情報も照合 |
| {a,b}-inverse-overflow / {a,b}-end-overflow | \(D=2^{-1074},e=0\)／\(D=2^{-1023},e=2^{-1024},\varepsilon_{\rm rel}=1/2\)。後者の中心は有限 | C `weight_length_not_representable / multiplier_not_representable`、α/β、stage=0。逆数と上端を区別 |
| {a,b}-release / {a,b}-ties | 保留からちょうど上限で解決／局所の同率を破って分割（下記） | 成功、分母・分割パネル・段数・呼出し履歴。保留中の畳み込み0回 |
| {a,b}-negative-hold | \(D=-2^{-12},e=|D|/2\)、初期1枚、上限2、相対条件を満たさない。次段で \(D=1,e=0\) | 成功。最初は負の分母を保留し、分割1回、次段で初めて確定。符号だけを見た入力拒否を検出 |
| {a,b}-up{0,1} | `_run`、該当 phase/stage **だけ** M3a の実例外を注入 | C `kernel_not_converged`、指定 phase と整数 stage。`inner_exception, inner_reason, inner_attributes, index, loss_eV, column, node_kind` が注入と一致 |
| c-dominates / exact-sum | αの \(A=2^{-40},e=2^{-39},c=1\)／\(A=1,c=3\cdot2^{-54},e=0\) | 成功。初段の13aの生予算は §4.4、分母条件は \(A+c\) で通る。最近接と外向き両端を §4.5 でビット比較 |
| other-point / other-loss | βの S2 分母確定、別列の点の誤差／別列の検証誤差だけが上限で予算超過 | C `kernel_error_exceeds_budget / loss_budget_not_met`、β、stage=null。初段の診断 `diagnostics.stages=0`、最悪の列と内訳を照合 |
| width / cells | βの次のメッシュで節点の分離不能／セル上限超過 | C `panel_width_unresolved / cell_limit_reached`、β、stage=null。1回二分したメッシュの評価前に失敗し、`diagnostics.stages=0`（最後に評価した初段）。パネル両端と幅／セル数と上限を照合 |

V9-pre の3行はすべて、公開 API で `weight_length_kind="inverse_13a_window_area", weight_length_A=None` を指定し、段αの開始回数（種の積分への進入を数える `seed_entries`）が0であること、`phase=null, stage=null`、核の呼出し・送信点数0を検査する。pre-N35a と pre-N32a の末尾はレビュー前の行 ID である。入力はこの表と JSON の `inputs` で指定したものを使う。

**段を選んで注入する方法。** αでは \(c=0\) の導出、βでは外部重みで S2 を基本とする。上流の片段失敗には導出＋S2 を使い、非対象 phase は正常に通す。能動の合成核で実際の核呼出しを通した後、面積集計の指定列とその局所誤差・5項だけを下記の値に置く。未指定の列はその段の予算内、`c-dominates` は通常の面積床 \(10^{-10}\) のまま判定する。分母保留用の小面積の試験では `area_epsabs=1`。成功段で \(A_{13a}^\beta=A_{13a}^\alpha\) とし、αの成功を調べる行は種入口で終了してよい。\(A\) の生予算と \(D\) の相対条件の両方を診断から読む。分母の条件を満たさない段で符号を確定しない。負値を拒否するには相対条件と \(D+e<0\) の両方を要する。

- `positive-limit` と `release` は初期1枚、\(D=2^{-12},e=D/2\)、`epsrel=1e-6`。生予算は通るが分母を保留する。前者は `panel_limit=1`、後者は2。後者の次段は \(D=1,e=0\) で上限ちょうどで成功。`zero-limit` は初期1枚、上限2、両段で \(D=e=0\)。最悪の比が元来 null の場合は null を保つ。通常の最悪キーがある場合はキー・予算・推定・比を保存する。
- `ties` は初期端 [2,3,5,7]、上限4、\(D=e=0\)、局所誤差はすべて0。最長で左の [3,5] を二分し、次段 \(D=1,e=0\) で成功。局所誤差優先の対照も同じ試験入口で、初段 \(D=2^{-12},e=D/2\)、分母の局所誤差を \((e,0,0)\) とし、短い [2,3] が先に選ばれることを確かめる。
- `up0/up1` は N29 の注入にならい、対象バッチ内の最初の内部点で `DepthAverageConvergenceError("budget_not_met", …)` を送出する。他の段には注入しない。stage 1 が必要な場合は対象 phase の stage 0 の分母だけを `release` と同じ条件で保留させる。両 phase の上限を独立に数える。段αの核失敗時にβの呼出しは0。最初のバッチ失敗時は `error_budget.outputs={}`、`max_ratio=null`。
- `budget_point/loss` は1枚・上限1、分母 \(D=1,e=0\)、別列 `term_13b` の面積0、局所 `kernel_point_error` または `verification` に \(2\,\text{area_epsabs}\) を与える。`refined_width` は \(W=[2,2+384\,\mathrm{ulp}(2)]\)、1枚から強制分割。子の幅192 ulpで内部節点が端に重なる。`refined_cells` は [2,4]、δ \(p=10\)、出力端 [0,20]、初期1枚で3セル、`max_cells=3`。これら2行は `release` の分母注入で第1細分を要求し、次段の評価前に幾何の失敗を観測する。

## 7. V10 の出力、照合、再生

`V10-{n}-{cli_area,cli_spectrum}` は short の API と CLI を同じ実効入力で走らせ、出口まで照合する。照合には、M4a の C02 と同じ形の小さい試験台帳を使う。入口ごとに short の**旧来の完全な conditions** を1件登録し、その入口の全キーに `pass` を置く。`partial` の対照は一次線または端だけを変え、`absent` は元台帳の運動量だけを変えたレコードを使う。これは照合表示の試験用で、物理参照の件数に入れない。API/CLI の両方に同じ台帳を渡す。

| \(n\) | 重み、\(c\)、乗数 | 模型 ID、群の照合 | 追加の警告と区間 |
| --- | --- | --- | --- |
| 1 | 外部20、乗数省略 | 旧、完全一致なら matched | 追加なし、両区間 null |
| 2 | 外部20、外部乗数20（独立に指定） | 旧、1 と同じ | 追加なし、両区間 null、乗数の出どころを保存 |
| 3 | 外部20、外部乗数7 | 新、not_matched | W-len、W-sens、両区間 null |
| 4 | 外部20、S2 | 新、not_matched | W-S2、W-sens、値が異なれば W-len。乗数のみ全量の逆区間 |
| 5 | 窓13a、\(c=0\)、同長 | 新、not_matched | W-13a、W-sens、両区間がαの逆区間 |
| 6 | 窓13a、\(c=1/32\)、外部乗数7 | 新、not_matched | W-13a、W-c、W-sens、異なれば W-len。重みのみαの逆区間 |
| 7 | 窓13a、\(c=1/32\)、S2 | 新、not_matched | W-13a、W-c、W-S2、W-sens、異なれば W-len。重みはα、乗数はβの逆区間 |

旧 ID は `simonsen1997_single_loss_spectrum_external_v0`、新 ID は `simonsen1997_single_loss_spectrum_conventions_v0`。各行の API の `model_id`、CLI の stdout、`run.json`、結果欄で同一の文字列を要求する。警告の W 記号は設計 §4.3 の**意味**を識別する。W-len の等値判定は報告した binary64 の長さと乗数で行う。3–7 行は `comparison.key_status` の各キーが `no_reference`、各群 `not_matched`、`physical_validity="not_established"`。1–2 行は旧版の警告列・照合結果そのものを要求する。`V10-{1,2}-{partial,absent}` は順に `kernel_conditions_checked, not_matched`。§5 の `V10-coincidence` は長さと乗数が一致しても行4のまま W-S2 は存在し、W-len は出ない。候補の照合条件は3–7で `weight_length_kind` と乗数の規約・値、導出時は \(c\) を含める。

- `run.json` の `spectrum_definition.weight` は `length_A,length_kind,length_source`。導出時に限り種、\(A^\alpha,e^\alpha,c\)、\(c\) の出どころ、分母条件、区間、α・βの整合の診断を持つ。`multiplier.{convention,value_A,interval_A,source}` と結果の対応欄を照合。行1–4の重み、行2/3/6の外部乗数、行6/7の \(c\) の出どころには互いに異なる入力文字列を渡し、どの欄に入ったかを比較する。導出種を省略した行5/7では20、行6では明示10 Å。長さ・端は Å、面積・誤差・\(c\) は Å⁻¹、感度の単位は §1。段αの設定は面積6キーのみ、上位の設定値・出どころは従来の6/8キー。入力、ハッシュ、模型・依存版、CSV、時刻、警告も V1 の欄を突き合わせる。
- 新形式の再生は、各 CLI の `run.json` だけを入力に、材料・初期エネルギー・窓・運動量・法線角・電荷・設定、一次線と出力端を復元する。導出時は `weight_length_A=None` と記録した種・\(c\)、S2 は S2 の指定を再投入する。保存した結果長や乗数を外部値で注入しない。復元した設定ファイルと（スペクトルの場合は）一次線・端のファイルを新しい一時先に書き、該当する CLI を再実行する。元の入力ファイルにはアクセスしない。双方の CSV バイト列、API の値・誤差・規約・区間・診断と対応する出どころを照合する。復元ファイルを再直列化した場合、その入力ハッシュは再構成した実バイトから別途検算する。API が受け取らない CLI 固有の出どころは記録同士で照合する。
- `V10-legacy-{cli_area,cli_spectrum}` は V1 の **基準版で生成した**旧 manifest を単独で読む。sdist では基準版の既存スキーマだけを持つ同じ入力の fixture を使える。新キーなしから `external_constant, multiplier=weight_length`、旧 `length_A` を補い、旧欄を再現。`V10-pos-{型名}` は基準の `dataclasses.fields` 順で親・子それぞれに旧位置引数だけを渡す。旧欄と子の並び、新7欄（旧 ID、外部重み、重み区間 null、\(c=0\)、同長規約、`multiplier_A==weight_length_A`、乗数区間 null）を検査する。
- `V10-caps` は `capabilities --json` の `xps_single_loss.model_id` が旧 ID で、同ブロックが両 ID・2種類の重み・3種類の乗数を記載すること、`supported` の全要素と順・その他のブロックが基準と一致することを確かめる。
- `V10-{cli_area,cli_spectrum}-{flag_case}` は次表の旗だけを有効な argv から変える。全行で終了2、stdout は空、stderr は単一の JSON（`status="error"`、`error_type="ValueError"` と文言）、出力先なし、核の呼出し0を要求する。`reason, phase, stage` は使い方の例外では該当なし。`V10-stderr-{a,b}` は導出＋S2、新規約を通る CLI に §6 の上流失敗を注入する。終了2、出力先なし、`status="not_converged", reason="kernel_not_converged"`、該当 phase、整数0、元の上流属性を照合する。

| flag_case | 誤った argv |
| --- | --- |
| no_weight / both_weights | 重みの旗を両方省略／`--weight-length-A 20` と `--weight-length-inverse-13a-window-area` を併記 |
| derived_source | 窓13aと `--weight-length-source text` |
| external_c_zero / external_c_source / external_seed | 外部重みに順に `--extra-inverse-length-per-A 0`、`--extra-inverse-length-source text`、`--weight-length-seed-A 20` |
| missing_c_with_source | 窓13aで \(c\) の旗を省略、`--extra-inverse-length-source text` |
| both_m / missing_m_with_source / s2_with_source | 乗数の両旗を併記／乗数省略で `--multiplier-source text`／S2 と同じ出どころの旗を併記 |

## 8. V11 の必須変異

全26件。設計 §6 の候補はすべて採用し、異なる出力が検出する誤りは別 ID にした。正常対照と1か所だけ誤らせた版を同じ条件で走らせる。指定キーで正常対照が合格、変異が不合格なら検出。起動失敗、他の行の失敗、参照不足、未確定、未実行に検出点を与えない。変異は全件プロファイルで対象行だけ再実行すればよい。

| V11 ID | 誤らせる操作 | 検出する行、キーと理由 |
| --- | --- | --- |
| M01 | 生の面積にも外部乗数を掛ける | V2-area-only、`window_area.total` のビットが変わる |
| M02 | S2 の分母を13aにする | V4-RB1、`multiplier_A` が最終全量の最近接逆数と異なる |
| M03 | 導出長に全量を使う | V5-RB1a-z-20、`weight_length_A` が13aの逆数の検算／参照区間から外れる |
| M04 | S2 の分母に \(c\) を足す | V2-api_spectrum-inverse_window_area、`cutoff_sensitivity.total`。同じ核の対照で乗数が変わる |
| M05 / M06 | \(c\) を中心／区間の端から落とす | V9-c-dominates、`weight_length_A / weight_length_interval_A`。\(c\) が支配する分母を有理数と照合 |
| M07 / M08 / M09 | \(H,w,\mathrm{LbK}_b\) にそれぞれ \(L_w\) を使う | V7-ext、`cutoff_sensitivity.total, single_loss_weight, diagnostics.LbK_b.total`。\(m=4,L_w=8\) |
| M10 | 収支の目標に \(L_w\) を使う | V7-ext、`diagnostics.balance.single_loss_total` の差が解析目標と合わない |
| M11 / M12 | 種の面積を返す／段βの核へ種を渡す | V6-wire、`window_area.term_13b`／`beta.call_kwargs.weight_length_A`。13b が長さで変化 |
| M13 | α・βの整合検査を外す | V6-outside、`exception.reason`。境の次の binary64 で成功してしまう |
| M14 | 逆数区間を内向きに丸める | V9-zero-c-weight_length、`weight_length_interval_A`。\(1/3\) を挟む隣接数が異なる |
| M15 | 最初に解決した S2 乗数を保持する | V4-stages、`retained_panel.folded_values`。未分割のセルが古い乗数を持つ |
| M16 / M17 | 衛星誤差の補正を省く／一次線にも加える | V4-error、`signed_error / zero_loss_error`。符号つき・非零幅の厳密例 |
| M18 | 分母の保留中に成功または停止する | V9-b-ties、`split_panel`。通常の分割集合が空でも最長・左を分割する |
| M19 | 全スペクトルを一律に2倍する | V3-RB1、`bin[6].spectrum`。このビンの一次線は1 |
| M20 | 新規約にも旧台帳の matched を返す | V10-3-cli_spectrum、`comparison.spectrum`。試験台帳に新規約の記録はない |
| M21 | 外部同値時だけ正の衛星ビンを1 ulp ずらす | V1-api_spectrum-equal_external、`single_loss_total` のビット |
| M22 / M23 | 同値外部の stdout に新 ID／既定の警告列に W-sens を追加 | V10-2-cli_area、`stdout.model_id`／V1-api_area-omit、`warnings`。旧経路の契約 |
| M24 / M25 | 左／右の外を全量から残りの領域を引いた残差で埋める | V3-tails、`outside.low.single_loss_total / outside.high.single_loss_total`。真の尾は2と4 |
| M26 | 分散あり parallel の13aをバルク API で置き換える | V5-RB2a-z-20、`seed.area13`。同じ運動量規約の R-B を Å⁻¹ で直接比較 |

## 9. 時間、凍結、変更

- CI は旧72行の同じ実行を再利用し、増分を各代表 CI 環境の壁時計で **60秒** までとする。全465行、初期 CI 236行（旧72を含む）。新しい実材料の CI は Al、1000 eV、[1,50]、`k_total` の一条件。V3/V4 はスペクトル、V5 は種の面積を1回。実行中に同一入力・設定・注入・版・入口の結果を共有してよい。独立再実行を要求する対、再生、異なる phase・規約の実行を1回で代用しない。
- 継承した CI 単独と拡張した CI を同じジョブ・並列度で測り、差を記録する。過負荷のため順送りする場合は `ci_reduction_order` のとおり **V4-RB1 → V3-RB1** の順に1行ずつ全件へ移し、移すたび再計測する（236→235→234）。V5 の Al、V3-layer と V4 の合成例、V1・V2、エラーと CLI は維持。手元の2回の最大を使うと最初の実材料3行は約 \(2\times16.1+13.0=45.2\) 秒。両方を全件へ移した場合は、移した2行の ID、各回の計測値、残った CI 234行を実行記録に残す。追加分で独立な R-B 参照と直接比べる CI 行は V5-RB1a-z-20 だけになり、V3・V4 は合成例で確かめる。CI の実測は実装後に記録する。
- 全件は直列1ワーカー、通常の行は2時間、変異全26件の対照と対象行の再実行を別枠で2時間とする。全件の種比較は同一の段α観測を複数の参照 ID と照合できる。上限に達したら残りを `not_run` として記録し、必須行を削らず再登録する。結果の値を開いてから行を選んだり許容を広げたりしない。
- 凍結候補で生成器を一度実行し、別のハッシュ種で `--check` する。入力ハッシュ、行 ID の一意性、CI 順、条件付き参照を含む予定キー、変異の対照と検出先をレビューする。統括担当が利用者の了解を記録して本文と状態表を凍結し、その版のコミットと付録のハッシュを記録する。その後に別担当が写像表を作って別途ハッシュで固定し、実装・参照の比較を開く。
- **凍結後の変更には利用者の了解が必要。** M4a §11 と同じく、旧文書・旧 JSON・旧ハッシュは保存する。変更台帳に「行・キー／凍結した期待／新しい期待／観測した環境と証拠／原因／同意日／影響する試験／再検証／新旧ハッシュ」を記す。必要な箇所だけ新版または明示した上書きで別試験として登録する。未確定への追試では参照担当へ必要なキーと精度のみを伝え、事前の精密化手順・旧値・新値を保存する。設計 U1–U5 の利用者決定に加え、統括担当が今回指定した例外属性・診断・誤差の丸めを本案に反映した。本文・表の正式凍結は別途行う。

### 設計の細部をここで決めた事項

| 設計第3版の箇所 | 事前登録で固定する細部 |
| --- | --- |
| §3.3・§4.2、§3.6 の互換経路 | §6 の `phase`・`stage` の割当て、継承する例外の属性、整合失敗の最終合格段。1–2行の `diagnostics` は基準と同じキー集合とする（§3.1） |
| §3.2 の構造的ゼロと負の分母 | 導出で \(c=0\) は重みを先に拒否する。\(c>0\) で長さが表せても S2 の全量は0なので `multiplier_undefined`。反復中の負の分母は \(D+e<0\) と相対条件の両方で確定する（§6） |
| §3.5、§6 V2・V4 | 窓13aの重みと乗数 `weight_length` の対は固定時の \(e_f\) を直接比べ、導出側の返却誤差を別に検算する（§3.2）。補正対象のキー、厳密な式の上向き丸め、密度の最近接丸めと S2 の \(w\) の丸めを §4.4 で固定した |
| §6 V3–V5・V8 と実行時間 | 固定乗数の写像は条件式を先に凍結する（§4.1）。派生密度と重複する13a 面積の集計は §4.4–4.5、V8 の再現幅は §5、CI の上限・削減順と記録は §9 に従う |

## 付録　凍結候補の識別

| 対象 | SHA-256 |
| --- | --- |
| M4a 状態表（入力、ファイル全体） | `250b45f863ccbd1c4c94474f23769976380c5dad5c2671328712ef4acdc7f408` |
| M4a R-B/CE（入力、ファイル全体） | `942353c748b4959220ba4da15f078a8ebdbdd9833d3abdc2e80cb8acfd6e8d94` |
| M4a-2 状態表（UTF-8、LF、末尾 LF を含むファイル全体） | `4a8264a58b9902efee3d843c10be812569ff6f46d720c9f23ebf37a4f263073f` |
| expected_states_m4a2.py | `b2feff2ea4b1a0efc097ec6a4103e5904531d8313081c6a10c57f183b24b11f2` |
| 本文（下記の自己参照欄正規化） | `4eccf1a10f27da8c6bbf0701be7e6b48cabe11581cedd050b2a8a75e445c7dea` |

本文のハッシュは、上表の「本文」行の64桁だけを ASCII の `0` 64個に置き換えたファイルバイト列（UTF-8、LF、末尾 LF あり）に対して計算する。ハッシュを書き込んだ後の**無加工の本文ファイル全体**の SHA-256 は別置きの報告にも記す。正式凍結時には同じ規則で再計算し、無加工のファイルのハッシュと版を外部の凍結記録に残す。
