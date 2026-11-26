# M4a-2 の設計　重みの長さと式(19)の乗数の選択肢

作成日: 2026-11-26。状態: **決定（利用者、2026-11-26、§9 の U1–U5 を推奨どおり）。第3版。事前登録の前。**計算のコードの基準は `e4e024f`（main `cf69b90` と `src/`・`tests/` が同じ）。`work/` の下のレビューと試算の記録は配布物に含めない。

## 作り方

- 出発点は M3a-2 の判定（[ROADMAP](ROADMAP.md) の「M3a-2 の成立性の判定」。作業の記録は `work/consult/m3a2/`）。利用者は 2026-11-26 に M4a-2 を選んだ。
- 統括担当（Claude Opus）が第1版を書いた。Codex（`work/consult/m4a2/codex_design_review.md`、blocker 3、major 6、minor 2）と Claude Fable（`fable_design_review.md`、blocker 0、major 8、minor 10）が互いを見ずにレビューし、どちらも「修正して採用」だった。第2版の確認で、Fable（`fable_design_recheck.md`）は「事前登録に進んでよい」（minor 8）、Codex（`codex_design_recheck.md`）は「直してから進む」（blocker 3、major 3、minor 4）だった。第3版はその全部を取り入れた（§11、§12）。
- 表記: [確認] コード・文書・文献で確かめたこと（担当を添える）、[導出]、[提案] 統括担当の案。レビューの試算の数値は設計の判断の材料で、検証の期待値には使わない（§6 の V8 の値は、M3a-2 の計算の記録 `work/consult/m3a2/compute/result_C-{Al,Si}.json` の値）。

## 1. 目的と範囲

M4a は、式(18)の深さの重みの長さと、式(19)の乗数に、同じ外部の長さ \(L\) を使う [確認: `xps_single_loss.py:196`–206 で `weight_length_A` を核の `kw` と `_loss_panels._integrate` の `length` の両方に渡す。`_integrate` の `length` は乗数としてだけ使われる（`_loss_panels.py:308`–318、350–354、463）。両担当]。M3a-2 の判定で、衛星の相対面積 \(w_{1,W}\) を大きく動かすのは、重みに使う一定長の選び方（Al 1000 eV で −11.9 %）と、式(19)の乗数の規約（Al で 1 から 1.52。1.52 は採らない M-P を含む）だと分かった。

M4a-2 は、この2つを呼び出す側が明示して選べるようにする。**物理は足さない。** 核（E6 の深さ平均）、窓の積分、セルの積分、誤差の予算は M4a のまま使う。

**範囲に入れる:** 重みの長さ \(L_w\) と乗数 \(m\) を別々に与えること。重みの長さの出どころに、外部の値（今と同じ）と、窓の 13a の面積から作る \(1/(A_{13a}+c)\)。乗数の規約に、重みの長さと同じ（今の既定）、外部の値、窓の面積の逆数（S2）。既定の呼び出しの旧来の数値は、M4a とビット単位で同じ。

**範囲に入れない:** 軌道ごとの重み（M3a-2 で本実装しないと決めた）、ポアソンの数え方の乗数 M-P（文献の式ではない。文書の注記だけ）、長さを決める窓と観測する窓を別にすること、\(\lambda_c\) の式の組み込み、\(c\) を足した乗数 \(1/(A_{\rm total}+c)\)（必要なら外部の乗数で与えられる。`window_area` は乗数によらない）、窓と規約を変えたときの全微分の感度、多重散乱、斜出射、窓の外の損失。

## 2. 定義

### 2.1 式 [導出。M4a の E7 で \(L\) を2つに分けただけ]

\[
J_{1,j}(E)=m\int_W F(E+T)\,K_{{\rm sc},j}(E_0,T;L_w)\,dT,\qquad
w_{1,W}=m\,A_{\rm total}(L_w),\qquad
H_j=m\,a\,K_{{\rm sc},j}(E_0,a;L_w)
\]

- \(A_j(L_w)\)（`window_area`、Å⁻¹）と \(J_0\) は乗数によらない。収支は \(\sum_DJ_{1,j}=S\,m\,A_j\)、\(\sum_DJ=S(1+mA_{\rm total})\)。
- **感度の意味** [Codex B1、確認の N1]: \(H_j=-\partial(mA_j)/\partial\ln a\)（無次元）、\(H^\pm_{j,D}=m\,a\,K_{{\rm sc},j}(a\pm0)\,G_D(a\pm0)=-\partial_\pm J_{1,j,D}/\partial\ln a\)（一次線の強度の単位）。どちらも**報告した \(L_w\) と \(m\) を固定した**偏感度で、窓を変えたときの \(m\) や \(L_w\) の変化は含まない（S2 では \(w\equiv1\) なので全微分は 0 だが、\(H\) は 0 でない。窓の長さでは \(dL_w/da=L_w^2K_{13a}(a)\) も効く）。窓を変えたときの値が要るなら、その窓で計算し直す。上端の診断 `diagnostics.LbK_b` は \(m\,b\,K_j(b)\)（無次元。キーの名前は変えない）。
- M4a の導出（`docs/M4A_DESIGN.md` §2.2）で \(L_w=m\) になるのは、一様な生成と一定長 \(L\) の指数の減衰を仮定したときだけである。\(L_w\neq m\) は、その導出から外れる規約の組み合わせで、導出の結果ではない。
- 一次線の強度が 0 のとき、S2 の \(w=1\) は規約の値で、衛星と一次線の観測の比（0/0）は定義されない。

### 2.2 重みの長さ \(L_w\)

| 種類（`weight_length_kind`） | \(L_w\) | 出典・根拠 |
| --- | --- | --- |
| `external_constant`（既定。今と同じ） | 呼び出す側が与える正の値 | Simonsen1997 p.1615 式(18) |
| `inverse_13a_window_area` | \(1/D\)、\(D=A_{13a}(W)+c\) | 出典なし。Tougaard2022 p.4 式(7)の形（面積と \(\lambda_c^{-1}\) の和の逆数）を、13a の項と有限の窓 \(W\) に当てはめたもの。式(6)–(7)の面積は環境の全有効断面積を 0 から \(E_{\max}\) まで積分したもので、この当てはめも重みへの使用も文献にない。根拠は M3a-2 の判定（軌道ごとの重み S1-A との差が小さかった） |

- 名前に `bulk` を使わない。`k_parallel` の 13a は一般に無限媒質のバルクの値と違う。どちらの運動量の規約でも「その規約の 13a 列の窓の面積」と定義する。
- 13a の被積分関数に深さも \(L\) も入らない [確認: `_depth_kernel.py:34`–40、`xps_depth_averaged.py:134`–173。Codex]。Al・Si、`k_total`・`k_parallel`、\(L=10\)、20、40 Å で `term_13a` の値と誤差はビット単位で同じ。窓の面積は \(L=5\)–80 Å で差が 2.1e-17 Å⁻¹ 以下 [確認: Fable の実測]。
- \(c\ge0\)（`extra_inverse_length_per_A`、Å⁻¹）は呼び出す側が与える追加の減衰。既定は 0。例えば \(\lambda_c^{-1}\)（Tougaard2022 式(8)）を入れられるが、パッケージは \(\lambda_c\) を計算しない（式(8)の係数の単位がこの頁だけでは確定できない）。\(c>0\) は M3a-2 で参考の条件（\(\lambda_c\)）でしか確かめていない。
- 意味は「同じ窓 \(W\) と同じ \(c\) のもとで、軌道ごとの重み S1-A に近い結果を与える一定長」に限る。M3a-2 の列挙した条件（Al と Si、1000 eV の主条件と参考の条件、法線出射、正孔電荷 1）で、面積と形の差が \(\Delta+u\le0.21\) % だった。重みの関数そのものの一致や、ほかの材料・電荷での精度は確立していない。**平均自由行程ではない。** 窓に依存する（Al で上端 50 eV と 99.8 eV で 27.06 Å と 25.96 Å）。外部から文献の IMFP を入れる使い方を否定しない。
- 長さを決める窓は、観測する窓 \(W\) と同じものに限る。外部の値の長さの出どころは M3a の呼び出しの境（`weight_length_A`）で、窓の長さでは M4a-2 がその値を作って M3a に渡す。

### 2.3 乗数 \(m\)

| 規約（`multiplier`） | \(m\) | 出典 | \(w_{1,W}\) |
| --- | --- | --- | --- |
| `weight_length`（既定。今と同じ） | \(L_w\) | Simonsen1997 p.1615 式(18)–(19)（同じ \(\lambda\)） | \(L_wA_{\rm total}\) |
| `external` | 呼び出す側が与える正の値 `multiplier_A`（Å） | 出典なし（比較のための規約の組み合わせ） | \(mA_{\rm total}\) |
| `inverse_window_area`（S2） | \(1/D\)、\(D=A_{\rm total}(W;L_w)\) | Yubero2005 p.045414-3 式(9)の後、Tougaard2011 pp.14、22 の規約（深さ平均の後の断面積の面積で乗数を決める）を、有限の窓 \(W\) に限って読んだもの | 定義により 1 |

- 文献の面積は、計算した損失の全範囲の面積である（Tougaard2011 p.14 は “usually taken”）。Yubero2005 は式(8)で軌道ごとの面積から経路の長さを、式(9)の後で深さ平均の後の面積から乗数を決める。Tougaard2011 p.21 は前者の自己無撞着な選び方を推奨する。
- S2 の \(D\) は窓 \(W\) の面積で、\(c\) を含まない（M3a-2 の計算と同じ）。全成分に同じ \(m\) を使う。\(c\) は重みを通して \(A_{\rm total}(L_w)\) に間接に効く。
- S2 では \(w_{1,W}\equiv1\) で、衛星の相対面積は予測ではなく定義になる。窓の中の形と、ビンごとの配分だけが意味をもつ。M-P は実装しない（一定長では \(L_w\) に一致する）。

### 2.4 組み合わせの表 [行の番号で規則を固定する。Codex 確認 M1]

| 行 | 重み | 乗数 | \(m\) | 模型の識別子 | `multiplier_convention` | 照合（`comparison`） | 追加の警告 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 外部 | 省略・`weight_length` | \(L_w\) | `simonsen1997_single_loss_spectrum_external_v0`（今のまま） | `weight_length` | M4a のまま | なし |
| 2 | 外部 | `external`、`multiplier_A` が \(L_w\) と binary64 で等しい | \(L_w\) | 同上 | `external`（出どころも記録） | M4a のまま（照合の条件は変えない） | なし |
| 3 | 外部 | `external`、等しくない | 与えた値 | `simonsen1997_single_loss_spectrum_conventions_v0` | `external` | `not_matched` | W-len、W-sens |
| 4 | 外部 | `inverse_window_area` | \(1/A_{\rm total}\) | 同上 | `inverse_window_area` | `not_matched` | W-S2、W-sens、\(m\neq L_w\) なら W-len |
| 5 | 窓の 13a | `weight_length` | \(L_w\) | 同上 | `weight_length` | `not_matched` | W-13a、W-sens、\(c>0\) なら W-c |
| 6 | 窓の 13a | `external` | 与えた値 | 同上 | `external` | `not_matched` | W-13a、W-sens、\(m\neq L_w\) なら W-len、\(c>0\) なら W-c |
| 7 | 窓の 13a | `inverse_window_area` | \(1/A_{\rm total}\) | 同上 | `inverse_window_area` | `not_matched` | W-13a、W-S2、W-sens、\(m\neq L_w\) なら W-len、\(c>0\) なら W-c |

- **1–2 行:** 旧来の数値の経路、旧来の識別子、旧来の照合の条件、旧来の警告の列をそのまま使う。照合は元の条件しだいで M4a の3つの状態を取りうる。2行目は計算が1行目と同じなので、識別子と照合も同じにし、利用者が独立に与えたことは規約・値・出どころの記録に残す。
- **3–7 行:** 新しい識別子を使い、照合の条件に `multiplier`（規約と値）と `weight_length_kind`（と、窓の 13a では \(c\)）を足す。登録の記録はすべて1行目の形なので、今の台帳では各群が `not_matched`、未照合のキーは `no_reference`、`physical_validity="not_established"` のまま。窓の長さや S2 の値が、たまたま登録の長さと等しくなっても、1–2 行に昇格させない。M4a の `kernel_conditions_checked` の意味（核の条件だけが一致し、一次線か出力の端が違う。`docs/M4A_PREREG.md:354`）は広げない。
- 新しい識別子は、式(19)の形に別の規約の長さを入れたものという意味で `simonsen1997_` を残す。文書にそう書く。
- 警告の記号: W-13a、W-c、W-len、W-S2、W-sens は §4.3。W-len は実際の \(m\neq L_w\)（binary64）で選ぶ。S2 の警告は、S2 の値がたまたま \(L_w\) と等しくても出す。
- API の結果、CLI の成功の JSON、`run.json`、`capabilities` の識別子は、同じ行で同じ値にする（今の CLI の後の2つは定数を直接使っている。`_single_loss_cli.py:195`、224。結果の `model_id` を使うように直す）。

### 2.5 分母と区間

逆数に使う量を**分母** \(D\) と呼ぶ。段αでは \(D=A^\alpha_{13a}+c\)、S2 の段βでは \(D=A_{\rm total}\)。その推定誤差を \(e\)（段αは `term_13a` の誤差、S2 は `total` の誤差。\(c\) は誤差なし）とする。

| 欄 | 外部の重み | 窓の 13a の重み |
| --- | --- | --- |
| `weight_length_interval_A` | `None` | \([1/(D+e),\,1/(D-e)]\) を外向きに丸めたもの。構造的なゼロで \(c>0\) なら \(e=0\) として \(1/c\) を外向きに丸めた区間（両端は \(1/c\) が binary64 で表せないとき違う） |
| `multiplier_interval_A`（乗数 `weight_length`） | `None` | 重みの区間と同じ |
| `multiplier_interval_A`（`external`） | `None` | `None` |
| `multiplier_interval_A`（S2） | \(A_{\rm total}\) の逆数の区間（M4a の `inverse_window_area_interval_A` と同じ値） | 同じ。報告した \(L_w\) のもとでの区間で、\(L_w\) の不確かさを通した同時の区間ではない |

区間は面積の**推定**誤差を変換したもので、証明された上界ではない。\(t=e/D\) とすると、報告した値から上端までの相対の距離は \(t/(1-t)\)、下端までは \(t/(1+t)\) で、上側が大きい。§3.2 の規則で \(t\le\varepsilon_{\rm rel}\) を満たすが、`epsrel` に上限はない（M4a の設定の検査は 1 未満を受け付ける）ので、\(t/(1-t)\) が小さいのは \(\varepsilon_{\rm rel}\) が小さいときに限る。

## 3. 数値の方法

### 3.1 計算の前の検査 [Codex M2、確認 B1]

窓、材料、設定、一次線、出力の端、初期のパネル数とセル数（M4a の N30 に相当）、新しい引数のすべてを、段αの前に検査する。**入力だけで決まる拒否**は、どれも核を評価する前に起こる。分母の符号による拒否（§3.2 の5）は、核を評価した後に決まる別の場合である。

### 3.2 分母の規則 [Codex B2・確認 B1、Fable M6・確認 minor 2、3]

1. **構造的なゼロ:** 強度が正の振動子がない材料は、核を呼ぶ前に判別する。
   - 窓の 13a の重みで \(c=0\): `SingleLossInputError("weight_length_undefined")`。S2 を同時に選んでいても、こちらを先にする。
   - 窓の 13a の重みで \(c>0\): \(A=e=0\) として \(L_w=1/c\)（最も近い binary64）と区間（§2.5）を作り、段αを飛ばす（`diagnostics["weight_length"]` に `seed_skipped="structural_zero"`、段αの核の呼び出し 0、面積と誤差は 0、整合の確認は `"not_applicable"`）。乗数 `weight_length`・`external` で \(J=J_0\)、\(w=0\) の成功。\(1/c\) か区間の端が表せなければ `weight_length_not_representable`。
   - S2 を外部の重みで選んだ: `SingleLossInputError("multiplier_undefined")`。
   - 活動線のある材料で面積がたまたま 0 になることを、入力の誤りとは扱わない。
2. **生の面積の予算は変えない:** 11列の面積には M4a の予算 \(\max(\varepsilon_{{\rm abs},A},\varepsilon_{\rm rel}|A_j|)\) をそのまま課す。`window_area_budget` と `diagnostics.error_budget` の `budget`・`ratio` も M4a の式の値のまま。
3. **分母の条件を足す:** 逆数を取る段では、さらに \(D-e>0\) と \(e\le\varepsilon_{\rm rel}|D|\) を課す。比較は、入力の binary64 を厳密な有理数にして行う（\(c\) との和、差、閾値を含む）。条件の値（\(D\)、\(e\)、\(\varepsilon_{\rm rel}|D|\)、満たした段）は `diagnostics["denominator"]`（段αは `diagnostics["weight_length"]` の中）に記録する。
4. **分母が決まらない段:** 条件を満たさない段では逆数を保留し、畳み込みと成功を禁止する（S2 の段βでは、その段の畳み込みをしない）。通常の分割の集合が空でも細分を続ける。分割するパネルは、分母の列の局所の誤差が最大のもの。同じ値なら最も長いもの、さらに同じなら左のもの。局所の誤差がすべて 0（\(D=e=0\) の段）でも、この同率の規則で分割する。
5. **符号:** 分母の条件の \(e\le\varepsilon_{\rm rel}|D|\) を満たし、\(D<0\)（\(D+e<0\)）の段で、`SingleLossInputError`（`weight_length_undefined`、`multiplier_undefined`）を返す。核を評価した後の拒否で、V9 では入力だけの拒否と別の行にする。条件を満たさない段では符号を決めない。
6. **上限と理由の優先:** `panel_limit` は**必要な分割が実行できないとき**に適用する（ちょうど上限の枚数で条件を満たした段は成功できる）。そのとき、分母の条件を満たしていなければ `SingleLossConvergenceError`（`weight_length_unresolved`、`multiplier_unresolved`）を返し、M4a の最悪のキーの記録も添える。分母の条件を満たし、ほかの予算だけが通らないときは M4a の理由（`kernel_error_exceeds_budget`、`loss_budget_not_met`）のまま。逆数か区間の端が有限の binary64 で表せなければ `weight_length_not_representable`、`multiplier_not_representable`。M4a のほかの理由（核の未収束、パネル幅の不能、セルの上限など）はそのまま。

### 3.3 窓の 13a の重みの2段の計算 [両担当が支持]

1. **段α（種の計算）:** 種の長さ \(L_{\rm seed}\)（`weight_length_seed_A`、省略で 20 Å）で、面積の積分を行う。設定は面積の6つ（`epsrel`、`area_epsabs`、`panel_limit`、`kernel_*`）だけを使う。**判定の対象は `term_13a` だけ**（M4a の面積の予算と §3.2 の分母の条件）。ほかの10列の値は捨て、予算も見ない [Fable 確認 minor 7]。ただし核は全列を計算するので、核の未収束は種によって起こりうる。種を公開の引数にするのは、そのとき種を変えられるようにするため。種は解析では結果に入らない。
2. **長さの確定:** \(A^\alpha_{13a}\)、\(e^\alpha\)、\(c\) を、足す前にそれぞれ厳密な有理数にする。\(L_w=1/D\) を最も近い binary64 に、区間の端 \(1/(D+e)\)、\(1/(D-e)\) を外向きに丸める。
3. **段β（本計算）:** \(L_w\) で M4a の計算（面積だけ、またはスペクトル）を行う。11列すべてを M4a の予算で判定する。返す面積は段βの最後のパネルから採る。
4. **整合の確認:** 報告した binary64 を厳密な有理数として \(|A^\alpha_{13a}-A^\beta_{13a}|\le e^\alpha+e^\beta\) を確かめる（任意の余裕は足さない）。外れたら `SingleLossConvergenceError("weight_length_inconsistent")`。役割は配線の誤り（窓や運動量の取り違え、段の取り違え）の検出で、基準は丸めの差より8桁以上ゆるい。数値の質は §6 の V2 と V5 で押さえる。
5. **段の記録:** 例外に `phase`（`"weight_length_seed"`、`"production"`）を、整数の `stage` と別に持たせる。`panel_limit` などの資源の上限は段ごとに数える。`diagnostics["weight_length"]` に、種、\(c\)、段αの面積と誤差、分母の条件、区間、核の呼び出しの数、パネル数、継承した警告、整合の確認の値（両方の面積と誤差、差、許容）を入れる。最上位の `kernel_calls` は段βだけ（M4a のキーの意味を保つ）。

### 3.4 S2 の1段の計算 [両担当が支持]

- 段βの各段で、まず `_areas(panels)` で面積を求める [確認: `_loss_panels.py:456` は `_fold` の 463 より先]。分母の条件を満たした段で、\(m\) を厳密な \(1/A_{\rm total}\) の最も近い binary64 に決め、その \(m\) で全セル（分割しなかったパネルを含む）の値、誤差の5項、局所の配分、予算を計算し直す。`_fold` は毎段、全パネルを畳み直すので、古い \(m\) は残らない [確認: `_fold.py:171`–223。両担当]。再利用するのは同じ \(L_w\) の核の点の値とパネルだけ。
- 最後の段の同じ \(m\) を、`_finish` の \(w\)、\(H\)、`LbK_b`、両側の外、密度、収支の目標と丸めの幅に渡す。報告する \(m\) と `window_area.total` は同じパネルから出るので、\(w=mA_{\rm total}\) は丸めの範囲で 1 になる（Fable の実測で 4 条件とも 1.0 か 0.9999999999999999）。
- `_budget` は段の間の出力の差を使わないので、最後の報告値を固定した M4a の成功の条件を保てる。判定の対象は面積の11列と、各領域の \(J_0\) と衛星の3成分で、\(J\)、\(w\)、\(H\) は派生量（M4a と同じ）。

### 3.5 誤差の意味 [Fable M7・確認 minor 4、8、Codex 確認 B2]

- **成功の契約と `diagnostics.error_budget`** は M4a のまま、報告した \(L_w\) と \(\hat m\) のもとでの推定誤差 \(e_f\) で判定する（判定の対象のキーは \(e_f\le B\)）。
- **乗数に区間があるとき**（S2 と、窓の 13a の重みで乗数 `weight_length`。§2.5）、乗数に比例する出力（衛星の3成分、\(H\)、\(H^\pm\)、両側の外の衛星、`spectrum` の衛星の部分）の**返す誤差**を
  \[e_{\rm ret}=e_f+\rho\,(|\hat y|+e_f),\qquad \rho=\max(\hat m-m_-,\,m_+-\hat m)/\hat m\]
  とする（\(m_\pm\) は外向きの区間の端。\(\rho\) は厳密な有理数で求めて上向きに丸める。上側が大きい）。これは \([\hat y-e_f,\hat y+e_f]\times[m_-/\hat m,\,m_+/\hat m]\) の四隅を包む保守の式で、符号つきの成分や \(\hat y=0\) にも使える。独立な誤差は仮定しない。`spectrum` の返す誤差は、衛星の部分の \(e_{\rm ret}\) に \(J_0\) の誤差と足し算の丸めを加えたもの。密度は補正した後のビンの誤差から一度だけ作る。面積と \(J_0\) はそのまま。予算には入れない（細分の後に1回足すだけ）。
- **返す誤差と予算の関係:** この場合 `*_error` は `diagnostics.error_budget` の `estimate` と一致しない。\(\rho\) と、足す前の \(e_f\) を読めるようにする（\(\rho\) は `diagnostics["multiplier"]` に記録し、足す前の値は試験の入口で読めるようにする。\(\rho\) だけでは足す前の binary64 を復元できないため）。\(\rho\le\varepsilon_{\rm rel}/(1-\varepsilon_{\rm rel})\) 程度なので、返す誤差は、\(\varepsilon_{\rm rel}\) が小さいとき予算のおよそ2倍以内に収まる。
- **\(w\) の誤差:** S2 では、返した binary64 の \(|w-1|\) を厳密に求めて上向きに丸めた値（定義により 1）。窓の 13a の重みで乗数 `weight_length` のときは \(e_{\rm ret}\) の式を \(w\) にも使う。ほかの規約では M4a の式の \(L\) を \(m\) に置き換えた \(m\,e_A+\varepsilon|w|\)。
- **重みの長さを通した核の変化は足さない:** 窓の 13a の重みでは、\(L_w\) の不確かさが核 \(K(L_w)\) を通して出力を動かす分を `*_error` に入れない。その感度は一般に 1 より小さいとは言えない（例えば \(m=L_w\) の \(H_{13a}=L_waK_{13a}(a)\) の対数の感度はちょうど 1。M3a-2 で測った条件では、長さ ±10 % で面積が数 %）。区間（`weight_length_interval_A`）と文書で示す。上の乗数の分（\(\rho\)）は足すので、乗数 `weight_length` の比例の部分は返す誤差に入る。

### 3.6 既定の経路を変えないこと [両担当]

- §2.4 の 1–2 行は同じ内部の経路を通り、M4a と同じ演算の列で同じ binary64 を返す。
- 1–2 行の照合の条件、`effective_settings`（面積6キー、スペクトル8キー）と出どころ、警告の列は変えない。新しい引数を `effective_settings` と照合の条件に漏らさない（`_single_loss_spectrum` の `**numerics` は `settings` に入る。`xps_single_loss.py:231`）。
- `integrate_window` と `fold` の署名は変えない（層の試験が使う）。乗数の選択肢を層の入口で試すなら、既定値つきのキーワードを足すだけにする。

## 4. API・CLI・記録

### 4.1 API [両担当の案を採った]

```python
# 外部の重み（今と同じ）。乗数を重みと別にするなら multiplier="external" と multiplier_A を与える
loss_window_area(mat, E0, loss_window_eV=(1., 50.), weight_length_A=19.4, epsilon_momentum="k_total",
                 multiplier="external", multiplier_A=24.6)
# 窓の 13a の重み。weight_length_A=None を明示する（必須のキーワードのまま）
single_loss_spectrum(mat, E0, primary, edges, loss_window_eV=(1., 50.), weight_length_A=None,
                     weight_length_kind="inverse_13a_window_area", extra_inverse_length_per_A=0.,
                     epsilon_momentum="k_total", multiplier="inverse_window_area")
```

署名に足すキーワード: `weight_length_kind="external_constant"`、`extra_inverse_length_per_A=0.`、`weight_length_seed_A=None`、`multiplier="weight_length"`、`multiplier_A=None`。

- 単位つきの数値の引数は数値だけにし、規約は文字列の引数で選ぶ。
- **検査の順** [Codex 確認 B3、Fable 確認 minor 1]: `_run` の今の長さの検査の位置で、まず `weight_length_kind` を確かめる（`str` か、列挙の値か）。`external_constant` なら、今の型・文・順の検査をそのまま行う（真偽値・NaN・無限大・0 以下は素の `ValueError`。`model.py:80`–86）。`inverse_13a_window_area` なら `weight_length_A is None` を確かめる。正の値の検査は確定した \(L_w\) に行う。ほかの新しい引数（\(c\)、種、乗数）の検査は、今の検査のあと。複数の引数が同時に不正なときの理由は、この順で決まる。
- `extra_inverse_length_per_A`: 0 以上の有限の実数。0 でない値は窓の 13a の重みのときだけ（API で明示した 0 は省略と同じ）。`weight_length_seed_A`: `None`（窓の 13a の重みでは 20 Å を補う）か正の有限の実数で、`None` 以外は窓の 13a の重みのときだけ。
- `multiplier`: `"weight_length"`、`"external"`、`"inverse_window_area"`。`multiplier_A`: `"external"` のときだけ正の有限の実数、ほかは `None`。§2.4 の2行目（等しい値）は `multiplier="external"` と `multiplier_A` で与える。
- 文字列の引数は、まず `str` か確かめてから列挙の値と比べる（配列や `np.bool_` との比較を避ける）。数値の文字列、真偽値、`np.bool_`、NaN、無限大、`OverflowError` を試験する。
- 結果に足す欄（`LossWindowAreaResult`。`dataclasses.field(kw_only=True, ...)` で既定値つきにし、子の型の位置引数の順を壊さない。Fable が Python 3.10 で確かめた）: `model_id`、`weight_length_kind`（既定 `"external_constant"`）、`weight_length_interval_A`（既定 `None`）、`extra_inverse_length_per_A`（既定 0）、`multiplier_convention`（既定 `"weight_length"`）、`multiplier_A`（既定 `None` で、`__post_init__` で `weight_length_A` を入れる）、`multiplier_interval_A`（既定 `None`）。`weight_length_A` は確定した \(L_w\)。旧来の位置引数だけで親・子の型を作っても整合する。
- 定数 `SINGLE_LOSS_MODEL_ID` は `_single_loss_types.py` に移し、`xps_single_loss.py` と `xpsloss` から同じ名前で再公開する（型の既定値に使うため。循環の読み込みを避ける）。値は変えない（`tests/m4a_states.py:1404`）。新しい識別子の定数を足す。

### 4.2 例外と理由

- 新しい `SingleLossInputError` の理由: `invalid_weight_length_kind`、`invalid_weight_length`（窓の 13a の重みで `weight_length_A` が `None` でない）、`invalid_extra_inverse_length`、`extra_inverse_length_without_derived_weight`、`invalid_weight_length_seed`、`weight_length_seed_without_derived_weight`、`invalid_multiplier`、`invalid_multiplier_length`、`weight_length_undefined`、`multiplier_undefined`。理由の文には、正しい引数の名前を入れる。
- 新しい `SingleLossConvergenceError` の理由: `weight_length_unresolved`、`weight_length_not_representable`、`weight_length_inconsistent`、`multiplier_unresolved`、`multiplier_not_representable`。
- `phase` を例外の属性に足す（§3.3）。`_ATTRIBUTES` に足すと、全例外の `vars()` と CLI の stderr の JSON にキーが増えるので、M4a の試験（`tests/m4a_states.py:913`、962）への影響を事前登録で確かめる。

### 4.3 警告

- **W-13a**（窓の 13a の重み）: 「この長さは窓に依存し、平均自由行程ではない」「軌道ごとの重みとの近さを確かめたのは Al と Si、法線出射、正孔電荷 1 の列挙した条件だけ」。
- **W-c**（\(c>0\)）: 「追加の減衰は呼び出す側が与えたもので、パッケージは計算していない。M3a-2 で確かめたのは参考の条件だけ」。
- **W-len**（実際の \(m\neq L_w\)）: 「重みの長さと乗数が違う。M4a の導出（一定長の減衰）から外れる規約の組み合わせ」。
- **W-S2**: 「この規約では単一損失の重みは定義により 1 で、衛星の相対面積の予測ではない」。
- **W-sens**（§2.4 の 3–7 行）: 「`cutoff_sensitivity` は報告した重みの長さと乗数を固定した偏感度」。
- M4a の警告はそのまま。1–2 行では警告の列を変えない。

### 4.4 CLI

| 旗 | 組 | 相手がないとき |
| --- | --- | --- |
| `--weight-length-A <値>` / `--weight-length-inverse-13a-window-area` | どちらか一方が必須 | 使い方の誤り（終了 2） |
| `--weight-length-source` | 外部の重みのときだけ | 同上 |
| `--extra-inverse-length-per-A`、`--extra-inverse-length-source`、`--weight-length-seed-A` | 窓の 13a の重みのときだけ（外部の重みでは値が 0 でも拒否）。出どころは \(c\) を与えたときだけ | 同上 |
| `--multiplier-A <値>` / `--multiplier-inverse-window-area` | どちらか一方、または省略（重みの長さと同じ）。`--multiplier-A` は API の `multiplier="external"` と `multiplier_A` に写す | — |
| `--multiplier-source` | `--multiplier-A` のときだけ | 使い方の誤り |

- 種は未指定で受け、窓の 13a の重みで 20 Å を補う（省略と無効な明示を見分ける）。
- `run.json` の `spectrum_definition.weight` は今の3キー（`length_A`、`length_kind`、`length_source`）を保つ。`length_kind` の値は `"external_constant"`（今と同じ）と `"inverse_13a_window_area"`。窓の 13a の重みのときだけ、13a の面積、誤差、分母の条件、区間、\(c\) とその出どころ、種、整合の確認を足す。`spectrum_definition.multiplier`（`convention`、`value_A`、`interval_A`、`source`）を足す。単位: 長さと区間は Å、\(A_{13a}\)、\(e\)、\(c\) は Å⁻¹。
- `run.json` だけから、選んだ規約ごとに計算を再生できるようにする。新しいキーのない旧来の `run.json` は、外部の重み・乗数 `weight_length` として読む [Codex 確認 N4]。
- `capabilities`: `xps_single_loss` のブロックだけを書き換える。`model_id` は今の識別子を残し、選べる規約と識別子を足す。`supported` の列には足さない（`tests/m4a_states.py:1056`）。
- 誤りのときの stderr の JSON、終了 2、出力先を作らないこと、は M4a のまま。

## 5. 入力の受理と拒否

§3.1、§3.2、§4.1、§4.2 のとおり。

## 6. 検証の骨子（軽い事前登録の素案）

物理を足さないので、検証は (1) 既定の経路が変わらないこと、(2) 同じ条件の別の呼び出しとの一致、(3) 独立な参照への写像、(4) 配線、(5) 拒否と未収束、(6) 記録、(7) 変異で組む。参照のファイル（R-A、R-B、状態の JSON）は変えない（凍結。読んで写像するだけ）。R-B と CE は上流の M3a を共有する。

| ID | 確かめること | 参照・基準 |
| --- | --- | --- |
| V1 既定の不変 | 同じ OS・依存の版・環境の中で、基準のコード（`e4e024f` の `src/xpsloss` を `git archive` で取り出し、子プロセスで走らせる。仕組みは `tests/m4a_states.py:634`–665。Git の履歴のある checkout で必須にし、展開した sdist での扱いは今の試験にそろえる）と変更後のコードを対にして走らせる。比べるのは基準のコミットにある欄とキーのすべて（`diagnostics`、警告の列を含む。binary64 は符号つきのゼロを含めビット列で）、CSV はバイト列、`run.json` は旧来のキーの値・型（時刻、経過時間、一時のパスは構造だけ）。新しい欄は「増えるだけ」を別に確かめる。§2.4 の 1–2 行（省略、`"weight_length"`、等しい `multiplier_A`）は、その引数のある入口（API の面積・スペクトル、CLI）だけで比べる。M4a の 100 行は今の終了状態と台帳で回す | 基準のコード |
| V2 同じ条件の別の呼び出し [Codex 確認 M2、Fable 確認 minor 5] | (a) 窓の 13a の重みの結果と、確定した \(L_w\) を外部の値で渡し、**同じ乗数の規約**（外部なら同じ値、S2 なら S2）と同じ入力・設定で呼んだ結果で、段βの数値、誤差、予算、パネルと呼び出しの数、収支などの決定的な数値の診断のキー（事前登録で列挙する）がビット単位で一致する。警告、照合、段αの診断、来歴、`model_id`、規約・区間の欄は除き、V10 で各実行の期待値と比べる。(b) 外部の乗数の面積だけの結果で、11列の面積、誤差、予算、そのキーの誤差の内訳が、乗数を省いた呼び出しとビット単位で一致する（`LbK_b`、\(H\)、\(w\) は乗数に従う）。V2 は配線の恒等式で、独立な照合は V3–V5 | 同じパッケージ |
| V3 外部の乗数の写像 [Codex B3] | R-B の \(L_w=L_{\rm ref}\) の行で \(q=m/L_{\rm ref}\) とし、キーごとに写像する。面積11列と \(J_0\) は参照のまま。衛星の3成分、\(H_j\)、\(H^\pm_{j,D}\)、\(w\) は \(q\) 倍。`spectrum` は \(J_0^{\rm ref}+qJ^{\rm ref}_{1,\rm total}\)（凍結した RB1 の `bin[6]` と RB3 の `outside.high` は衛星が 0 で、\(q=2\) でも 1 と約 0.35 のまま）。密度は写像した値をビン幅で割る。両側の外もビンと同じに直接比べる。許容は `docs/M4A_PREREG.md` §7.2 の式を \(L\to m\) で作り直す（\(r\to qr\)、\(u_r\to qu_r\)、\(d_{\rm rel}\) はそのまま、\(d_{\rm abs}\) はビンと外が \(\eta S_{\rm abs}\) のまま、\(w\) は \(m\) × 面積の床、\(H\) は \(m\,a\,10^{-12}\)、\(H^\pm\) は \(m\,a\,10^{-12}S_{\rm abs}\)）。和の不確かさは \(u_0+qu_1\) と丸め。判定は M4a の精度の基準 \(|y-r|+u_r\le d\) の三値 | R-B（独立な参照） |
| V4 S2 [Codex 確認 B2、M3] | ① 足す前の \(e_f\) を、報告した \(\hat m\) を固定した写像と三値で比べ、判定の対象のキーで \(e_f\le B\)。② 返す \(e_{\rm ret}\) を、独立な分子と分母の参照から区間の除算で作った S2 の比と三値で比べる（\(q=1/(L_{\rm ref}A^{\rm ref}_{\rm total})\)。分母の下端が正であることを確かめる。符号つきの `hole` は四隅で囲う。相関の独立は仮定しない。床には参照から得た \(m_{\rm ref}\) を使う）。③ 足し算と丸めの式そのもの（§3.5）を検算する。\(|w-1|\) は丸めの恒等式として数える。全領域の衛星の和が \(S\)。参照の不足と区間の部分的な重なりは別に数える | R-B |
| V5 窓の 13a の長さ [Fable M3、Codex 確認 M3] | R-B の `window_area.term_13a` の参照（18 行。RB 系 12、CE 系 6）で、まず段αの \(\hat A^\alpha_{13a}\) を Å⁻¹ で M4a の精度の基準 \(|\hat A-r|+u_r\le d\) と照合する。次に、報告した \(A\)、\(e\)、\(c\) の厳密な和から、最も近い \(L_w\) と外向きの両端を独立に検算する（\(c=0\) と \(c>0\)）。報告した区間の質は、参照の \(r\pm u_r\) を反転した区間との包含・分離・部分的な重なりを別に数える。種の依存（10、20、40 Å）は全件だけで、各種が参照と整合し、種どうしの区間が交わる | R-B |
| V6 配線 [Codex M6] | 13a が一定で、ほかの列が長さで変わる既知の合成の核を `_run` に通し、`_call_kernel` の実引数を観測する（段αで種、段βで \(L_w\)。\(L_{\rm seed}\neq L_w\) の条件で）。整合の確認に、許容の内側・等号・すぐ外の擾乱を入れる。実材料では、分散のある両運動量の規約を含め、`k_parallel` で 13a を別の量に取り替える変異を落とす | 合成の核 |
| V7 感度の意味 [Codex B1、確認 N1] | 正の定数の核と S2 で、2つの窓で \(w=1\)、各窓で \(H=a/(b-a)\) の 0 でない値。ビンでは \(G_D(a-0)\neq G_D(a+0)\) の境で左右の値を確かめる。同じ核で \(m\neq L_w\) として、`LbK_b` を \(m\,b\,K_j(b)\) と直接比べる | 解析の値 |
| V8 M3a-2 の再現 [Codex 確認 N3] | Al 1000 eV、窓 [1, 50] eV、`k_total`、\(c=0\): \(L_w\) が 27.0571259020694 Å（M3a-2 の `result_C-Al.json` の `B_length`）、乗数 19.4 Å で \(w=1.0861615888385654\)（同じファイルの `checks.B.M4a` × 19.4）。Si 窓 [3, 50] eV: 25.549842388377478 Å、乗数 24.6 Å で 1.6938314906967444。許容は、記録した入力と精度、M3a-2 の種（19.4、24.6 Å）と既定の 20 Å の違いから、凍結の前に決める。M3a-2 の表の B の値（1.086146、1.693802）は固定深さの数値の平均を含むので使わない。つなぎ方の確認で、独立な照合ではない（Si の [3, 50] eV は R-B にない） | M3a-2 の計算 |
| V9 拒否と未収束 | §3.2、§4.2 の各理由。入力だけの拒否は核の評価0回（N30 に相当）。分母の符号の拒否は別の行。段αだけ・段βだけの上流の失敗と `phase`・`stage`。構造的なゼロ（\(c=0\)、\(c>0\)、S2、\(c=0\) と S2 の同時）、強度0の一次線と活動のある材料、符号つきの成分、正だが決まらない分母（保留から細分に進むこと、\(D=e=0\) の段を含む）、ちょうど上限の枚数での成功、表せない逆数と区間の端、`panel_limit`、節点の分離の不能、セルの上限。\(c>0\) と S2 を同時に選んだ成功の例。API の検査の順（§4.1） | — |
| V10 記録 | §2.4 の7行それぞれの識別子（API、CLI の成功の JSON、`run.json` を同じ行で突き合わせる）、規約、警告、照合の状態。`run.json` の新しいキーと単位、両 CLI の再生（新しい形式と旧来の形式）、旗の組、stderr の JSON と終了 2。結果の型を旧来の位置引数で作ったとき、補った `multiplier_A` が `weight_length_A` と等しいこと | — |
| V11 変異 | 必須の ID を事前登録で固定し、各変異を落とす行・キー・理由を凍結する。候補: 乗数を `window_area` にも掛ける。S2 に \(A_{13a}\) を使う。窓の長さに \(A_{\rm total}\) を使う。\(c\) を S2 に足す。\(c\) を長さか区間に足し忘れる。\(H\)、\(w\)、`LbK_b`、収支の目標に \(L_w\) を使う。段αの面積を出力に使う。段βの核に確定した \(L_w\) の代わりに種を渡す（\(L_{\rm seed}\neq L_w\) の V6 で落とす）。整合の確認を外す。区間を内向きに丸める。S2 の \(m\) を最初に分母が決まった段で固定する（分母の決まった段が2回以上あり、\(m\) が変わり、分割しないパネルも残る行で落とす）。S2 の返す誤差への足し算を忘れる。\(J_0\) の誤差にも足す。分母の保留のときに細分せず終わる。スペクトル全体を \(q\) 倍する。3–7 行で `matched` を返す。2行目で既定と違う数値を返す。1–2 行で識別子か警告を変える。両側の外を残差で埋める（収支だけでは落ちない）。`k_parallel` の 13a をバルクの API で置き換える | 各変異を落とす行・キー・理由 |

- **時間** [Fable M8、Codex N2]: 実材料の面積の計算1回は、作業中の手元の機械で Al 3–10 秒、Si 5–18 秒だった（Fable の実測。CI の値ではない）。窓の 13a の重みでは、これがおよそ2倍になる（段αは 13a だけを判定するので、細分は少ないことがある）。CI と全件を分け、減らす順番を事前登録で決める。CI は Al の1条件と合成の核に絞り、乗数の写像は層の入口でも試す。種の依存は全件だけ。CI に足す時間の目安は、手元と CI で測ってから決める。

## 7. リスクと未確定の点

1. **\(L_w\neq m\) の物理の意味:** 導出の外の組み合わせ。規約を選ぶ道具として出し、どれが正しいかは主張しない。
2. **返す誤差:** 乗数に区間があるとき、予算の約2倍になりうる（\(\varepsilon_{\rm rel}\) が小さいとき）。`*_error` と `error_budget` の `estimate` が一致しない。文書に書く。
3. **\(c\) の扱い:** \(c>0\) は M3a-2 の参考の条件でしか確かめていない。
4. **13a の意味の違い:** `k_parallel` の 13a は無限媒質のバルクの値と違う。名前と文書で区別する。
5. **費用:** 窓の 13a の重みは呼び出しの時間がおよそ2倍。
6. **`phase` の属性:** M4a の例外の記録のキーが増える。

## 8. 選ばなかった方法

- **2回の呼び出しの手順を文書に書くだけ:** 乗数を重みと別にできない（今の API では \(m=L_w\) に固定）。
- **S2 の2段の計算（面積の段で \(m\) を決めてからスペクトル）:** 最後の \(m\) と報告する面積が別のパネルから出る。
- **区間の不確かさを予算に入れる:** ビンの予算に足すと予算を使い切りやすい。§3.2 で分母の条件を課して細分し、§3.5 で返す誤差に足す形を採った。
- **分母の列の予算を相対だけに置き換える（第2版）:** 生の面積の予算の記録が M4a と違う値になる。M4a の予算は残し、分母の条件を足す形に改めた（Codex 確認 B1、Fable 確認 minor 3）。
- **核の条件への射影で `kernel_conditions_checked` を出す（Codex 第1版 M3 の案）:** 凍結した語の意味を広げる。
- **`weight_length_A` と乗数の引数に数値と文字列を混ぜる（第1版）:** 単位つきの名前に文字列が入り、取り違えやすい。
- **段αでも11列すべてを判定する（第2版）:** 捨てる列の予算のために細分と失敗が増える。段αは 13a だけを判定する（Fable 確認 minor 7）。種の引数は、核の未収束のときのために残す。
- **\(\lambda_c\) の式の組み込み、長さを決める窓の分離、M-P、\(c\) を足した乗数:** §1 のとおり。

## 9. 利用者の決定（2026-11-26 に決定）

両担当とも、推奨に反対はなかった。利用者は 2026-11-26 に U1–U5 をすべて推奨どおりに決めた。あわせて D10 の出射角は「法線出射だけ」（M4a-2 の後は M4b を先にする候補）。

| ID | 決めたこと |
| --- | --- |
| U1 | 乗数に外部の値と S2 を入れる。M-P は入れない（文書の注記だけ） |
| U2 | 窓の 13a の長さを組み込みの選択肢にする（§3.3）。費用は呼び出しの時間がおよそ2倍 |
| U3 | 追加の減衰 \(c\) を外から与えられるようにする（既定 0、警告つき）。\(\lambda_c\) の式は組み込まない |
| U4 | 既定の振る舞いは M4a のまま。旧来の数値の欄はビット単位で保つ。結果と `run.json` の欄・キーは足すだけ |
| U5 | 軽い事前登録（§6。参照のファイルは新しく作らず、R-B の写像と、同じ条件の別の呼び出しとの一致で確かめる） |

## 10. 規模 [試算]

- 事前登録（Codex が下書き、統括担当と Fable がレビュー）: 0.5–1 人日。写像した参照の表（実装と別の担当）: 0.5 人日。
- 実装（Codex、worktree）: 1 人日。レビュー（Fable）と変異: 0.5–1 人日。文書: 0.25 人日。
- 合わせて約 3–4 人日。

## 11. 第1版のレビューの指摘と対応

| 指摘 | 対応 |
| --- | --- |
| Codex B1（感度の意味） | §2.1 に偏感度と明記。V7 |
| Codex B2（分母が決まらない段） | §3.2 の規則 |
| Codex B3（スペクトルの一律の倍） | V3 のキーごとの写像、V4 の区間の除算 |
| Codex M1（種、整合、区間の意味） | §3.3 の有理数の整合、§2.5、§3.2 の1 |
| Codex M2（計算の前の拒否、段の記録） | §3.1、§3.3 の `phase` |
| Codex M3（既定の互換、照合） | §2.4、§3.6、V1。照合は Fable M2 の案（3–7 行は一律に `not_matched`） |
| Codex M4（型、排他、結果の型） | §4.1、§4.4 |
| Codex M5（文献の規約） | §2.2、§2.3 |
| Codex M6（配線、V5 の許容） | V6、V8 |
| Codex N1（上端の診断、単位）、N2（費用） | §2.1、§4.4、§6 の時間 |
| Fable M1（V1 の SHA-256） | V1 を同じ環境の対の実行に |
| Fable M2（照合の規則） | §2.4 |
| Fable M3（V4 は R-A で実行できない） | V5 を R-B に |
| Fable M4（2回の呼び出しとの一致） | V2 |
| Fable M5（識別子と規約の規則） | §2.4、結果の `model_id` |
| Fable M6（絶対の床、負の面積） | §3.2 |
| Fable M7（S2 の返す誤差） | §3.5 |
| Fable M8（時間、V2 の許容） | §6 の時間、V3 の許容 |
| Fable minor 1–10 | §2.2、§2.3、§3.2、§3.3、§3.6、§4.1、§4.2、§4.4 |

## 12. 第2版の確認の指摘と対応

| 指摘 | 対応 |
| --- | --- |
| Codex 確認 B1（分母の定義、保留から進めること、構造的なゼロで \(c>0\)） | §2.5 で \(D=A_{13a}+c\)、\(D=A_{\rm total}\) を定義。§3.2 の3（分母の条件は有理数で \(D-e>0\)、\(e\le\varepsilon_{\rm rel}|D|\)）、4（保留中の分割の規則と同率の規則）、6（上限は必要な分割が実行できないとき）、1（\(c>0\) の構造的なゼロ）。§3.1 を「入力だけで決まる拒否」に限った |
| Codex 確認 B2（S2 の返す誤差と V4） | §3.5 の \(e_{\rm ret}=e_f+\rho(|\hat y|+e_f)\)、\(\rho\) は上側。V4 を ①②③ に分けた。「約2倍」は \(\varepsilon_{\rm rel}\) が小さいときに限った。感度が 1 より小さいという一般の主張を消した |
| Codex 確認 B3（`None` と検査の順） | §4.1 の検査の順と例 |
| Codex 確認 M1（§2.4 の表を共通の規則に） | §2.4 に行の番号、照合、警告の列。識別子は結果の `model_id` から出す |
| Codex 確認 M2（V2 の範囲） | V2 |
| Codex 確認 M3（V5 の参照の不確かさの向き） | V5、V3 の三値の判定 |
| Codex 確認 N1（感度の単位、V7 の境、`LbK_b`、0/0） | §2.1、V7 |
| Codex 確認 N2（V11 の変異の向き、足す変異） | V11 |
| Codex 確認 N3（V1・V8 の再現の条件） | V1、V8 の値（Codex が保存した値から確かめた） |
| Codex 確認 N4（旧来の形式からの再生） | §4.4、V10 |
| Fable 確認 minor 1（検査の順） | §4.1 |
| Fable 確認 minor 2（終わり方） | §3.1、§3.2 の1、5、6 |
| Fable 確認 minor 3（予算の記録） | §3.2 の2、3（M4a の予算を残し、分母の条件を別に記録） |
| Fable 確認 minor 4（\(\rho\) の定義、V4 のどちらの誤差か） | §3.5、V4 |
| Fable 確認 minor 5（V2 の範囲） | V2 |
| Fable 確認 minor 6（定数の場所、`None` の明示） | §4.1 |
| Fable 確認 minor 7（段αの判定の対象） | §3.3 の1（採った）、§8 |
| Fable 確認 minor 8（窓の重みと乗数 `weight_length` の誤差） | §3.5（乗数の区間の分を足す。核を通した分は足さない） |
