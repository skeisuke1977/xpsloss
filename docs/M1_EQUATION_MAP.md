# 式の対応表（Simonsen, Yubero, Tougaard 1997）

作成日: 2026-11-19。対象: A. Cohen Simonsen, F. Yubero, S. Tougaard, Phys. Rev. B **56**, 1612–1619 (1997)（論文庫キー `Simonsen1997`、[R2]）。ページ番号は誌面のページ（PDFの1–8頁が誌面1612–1619頁）。

## 1. 転記の方法と結果

- **方法:** 2つの担当（Claude Fable と Codex）が、互いの結果を見ずに、ページ画像（250 dpi）と400–450 dpiの切出しから転記した。PDFの抽出文字と二次資料は使っていない。
- **結果:** 番号付きの20式すべてについて、記号、符号、因子、べき、指数、括弧の範囲が2つの転記で一致した。判読不能の箇所はない。
- **照合:** 統括担当（Claude Opus）が、式(13)–(16)を切出し画像で直接照合した。
- **区別:** ここで確定したのは「何が印刷されているか」である。印刷された式が正しいか（誤植がないか）は、式レビューA（[REVIEW_CHECKLISTS.md](REVIEW_CHECKLISTS.md) のA1–A11）で別に検査した。
- **式レビューA（2026-11-19）:** 原典とこの表だけを渡された別の担当（Claude Opus。統括担当とは文脈を共有しない）が行った。式(13)–(16)は、θ = 0 で、式(2)、(3)、(6)、(8)–(12)と k⊥ の範囲の拡張から誤植なしに導かれた。k⊥ を数値積分する独立な経路と閉形式の相対差は、全項で 3.1×10⁻¹⁰ 以内だった。わざと1か所を変えた式は、すべてこの比較で検出された。報告は `work/reviews/M1_review_A_opus.md`、計算は `work/reviews/scratch_A/`（gitの管理外）。§5の Q4、Q5、Q14 はこのレビューで解決した。人による査読ではない。
- 転記の原本は `work/transcription/A_fable.md`、`A_codex.md` にある。gitの管理外で、公開しない。

速度は表示式で ν に似た字形で印刷されているが、本文の定義（p. 1613右、p. 1614右、p. 1615左）から速度 v とその成分である。以下では v と書く。誘電関数の字形は ϵ。

## 2. 原著の規約

| 項目 | 内容 | 位置 |
| --- | --- | --- |
| 単位系 | 明示されない。ポアソン方程式に 4π/k² が現れるガウス型の表記で、ħ、m、e、a0 は残っている | p. 1612–1615 |
| フーリエ変換 | 順変換 ρ(k,ω) = ∫dt∫d³r ρ(r,t) e^{−ik·r+iωt}、前因子なし。場は e^{ik·r−iωt} で戻す | 式(2)、(6) |
| 時間依存の規約 | 上の変換は e^{−iωt} の規約に当たり、式(20)は ω > 0 で Im ε > 0 を与える。本パッケージの e^{−iωt}、F = −1/ε の規約（[DIELECTRIC_DESIGN.md](DIELECTRIC_DESIGN.md)）と原著の符号の向きは一致する（式レビューA8で確認） | 式(6)、(20) |
| 座標 | 固体 r⊥ < 0、真空 r⊥ > 0、表面 r⊥ = 0。法線は固体から真空へ向く | p. 1612右–1613左、図1 |
| 深さ | 生成点 **a** = (−a, 0)、深さ a > 0。依頼書や manual の x0 は a に当たる | 式(3)の直後 |
| 角度 | θ は表面法線からの角度。法線出射が θ = 0 | p. 1613左 |
| 経路長 | 固体内の経路長 x = a/cos θ。K_eff はこの x 当たり | p. 1613左、p. 1615左 |
| 電荷と運動 | t = 0 に電子 −e と正孔 +e が同じ点に突然生じる。正孔は静止し寿命は無限。電子は等速直線で、t_s = a/v⊥ に表面を通る。等速の近似は ħω ≪ E0 で成り立つ。弾性散乱は無視する | 式(3)、p. 1613左、p. 1615左 |
| 境界 | 鏡面反射模型。疑似媒質と疑似真空、鏡像電荷、仮想表面電荷 σ を使う | p. 1613右、式(4)–(7) |

## 3. 記号と eV・Å での換算

H = ħ²/2m = 3.80998212 eV Å²、a0 = 0.529177 Å。長さは Å、エネルギーは eV、波数は Å⁻¹。

| 記号 | 意味 | eV・Åでの計算 |
| --- | --- | --- |
| E0 | 電子の初期運動エネルギー（非相対論、E0 = ½mv²） | 入力 |
| ħω | 損失エネルギー | 入力。0 < ħω < E0 |
| a | 生成深さ | 入力 |
| θ | 表面法線からの出射角 | 入力。M2-1では0のみ |
| x | 固体内の経路長 | a/cos θ |
| v⊥, v∥ | 速度の法線成分と面内成分 | v cos θ、v sin θ（導出による） |
| E⊥ | 法線方向の運動エネルギー ½mv⊥² | E0 cos²θ |
| ω⊥ | ω/v⊥（波数の次元） | ħω / (2 (H E0)^{1/2} cos θ)。ħv = 2(H E0)^{1/2} を使う |
| Ω0 | (ω − k∥v∥)/v⊥（波数の次元） | ω⊥ − k∥ tan θ |
| k± | 運動学的な全波数の上下限 | (√E0 ± √(E0 − ħω)) / √H |
| k∥± | 面内波数の積分範囲 | 式(17) |

## 4. 式と実装の対応

### 4.1 式（印刷のまま）

以下は転記を検証するために A. Cohen Simonsen, F. Yubero, S. Tougaard, Phys. Rev. B **56**, 1612 (1997) から組み直して引用した式であり、本プロジェクトの CC BY 4.0 の許諾対象外である。

式(1)（p. 1612左）と、直後の本文の定義:

$$
K(E_0,\hbar\omega)=\frac{1}{E_0\pi a_0}\int_{k_-}^{k_+}\frac{dk}{k}\,\mathrm{Im}\!\left[\frac{1}{\epsilon(k,\omega)}\right],\qquad
k_\pm=(2m/\hbar^2)^{1/2}\left[E_0^{1/2}\pm(E_0-\hbar\omega)^{1/2}\right].
$$

式(12)の後の本文（p. 1614）: d³k = 2π k∥ dk∥ dk⊥ とし、k⊥ の範囲を −∞ < k⊥ < ∞ に広げて k⊥ を解析的に積分する。その結果を4項の和として書く。

$$
K_{\rm eff}(E_0,\hbar\omega,a,\theta)=K_{M1}+K_{M2}+K_{M3}+K_V .
$$

式(13):

$$
\begin{aligned}
K_{M1}={}&\frac{\cos\theta}{\pi a_0E_\perp}\,\mathrm{Re}\!\left[i\int k_\parallel dk_\parallel\,\frac{(1/\epsilon-1)}{\Omega_0^2+k_\parallel^2}\right]
+\frac{-1}{2\pi E_\perp a_0x\omega_\perp}\,\mathrm{Re}\!\Bigg[i\int k_\parallel dk_\parallel\,\frac{(1/\epsilon-1)}{(\Omega_0^2+k_\parallel^2)^2}\\
&\times\Bigg\{\left(i+\frac{v_\parallel}{v_\perp}\right)(\Omega_0+ik_\parallel)^2\left(e^{-i\Omega_0a-k_\parallel a}-1\right)
+\left(-i+\frac{v_\parallel}{v_\perp}\right)(\Omega_0-ik_\parallel)^2\left(e^{i\Omega_0a-k_\parallel a}-1\right)\\
&\quad+\left(i-\frac{v_\parallel}{v_\perp}\right)\left(e^{-i\Omega_0a}-e^{-k_\parallel a}\right)\left(e^{i\Omega_0a}-e^{-k_\parallel a}\right)(\Omega_0^2+k_\parallel^2)\Bigg\}\Bigg]
\end{aligned}
$$

式(14):

$$
\begin{aligned}
K_{M2}={}&\frac{-1}{2\pi E_\perp a_0x\omega_\perp}\,\mathrm{Re}\!\Bigg[\left(\pi\delta(\omega_\perp)+\frac{i}{\omega_\perp}\right)\int k_\parallel dk_\parallel\,\frac{(1/\epsilon-1)}{\Omega_0^2+k_\parallel^2}\\
&\times\Bigg\{\left(i+\frac{v_\parallel}{v_\perp}\right)(\Omega_0+ik_\parallel)\left(1-e^{-k_\parallel a-i\Omega_0a}\right)
+\left(-i+\frac{v_\parallel}{v_\perp}\right)(\Omega_0-ik_\parallel)\left(e^{-2k_\parallel a}-e^{-k_\parallel a-i\Omega_0a}\right)\Bigg\}\Bigg]
\end{aligned}
$$

式(15):

$$
\begin{aligned}
K_{M3}={}&\frac{-1}{2\pi E_\perp a_0x\omega_\perp}\,\mathrm{Re}\!\Bigg[i\int k_\parallel dk_\parallel\,\frac{(-i+v_\parallel/v_\perp)(\Omega_0-ik_\parallel)\left(e^{-k_\parallel a}-e^{-i\Omega_0a}\right)}{\Omega_0^2+k_\parallel^2}\\
&\times\Bigg\{\frac{(\epsilon-1)(\epsilon+2)}{\epsilon(\epsilon+1)}\Bigg[\left(e^{i\Omega_0a}-e^{-k_\parallel a}\right)\frac{\Omega_0+ik_\parallel}{\Omega_0^2+k_\parallel^2}
-i\left(\pi\delta(\omega_\perp)+\frac{i}{\omega_\perp}\right)e^{-k_\parallel a}\Bigg]
+\frac{(\epsilon-1)}{(\epsilon+1)}\,\frac{\Omega_0-ik_\parallel}{\Omega_0^2+k_\parallel^2}\,e^{i\Omega_0a}\Bigg\}\Bigg]
\end{aligned}
$$

式(16):

$$
\begin{aligned}
K_V={}&\frac{1}{2\pi E_\perp a_0x\omega_\perp}\,\mathrm{Re}\!\Bigg[\int k_\parallel dk_\parallel\,\frac{(1-\epsilon)}{(1+\epsilon)}\,\frac{(1-iv_\parallel/v_\perp)(\Omega_0+ik_\parallel)e^{-i\Omega_0a}}{\Omega_0^2+k_\parallel^2}\\
&\times\Bigg\{\frac{2ik_\parallel e^{i\Omega_0a}-(\Omega_0+ik_\parallel)e^{-k_\parallel a}}{\Omega_0^2+k_\parallel^2}
-i\left(\pi\delta(\omega_\perp)+\frac{i}{\omega_\perp}\right)e^{-k_\parallel a}\Bigg\}\Bigg]
\end{aligned}
$$

式(16)の直後で E⊥ = ½mv⊥²、ω⊥ = ω/v⊥ と定義される。積分の上下限は k² = k∥² + Ω0² を k∥ について解いて得る（p. 1615左）。式(17):

$$
k_{\parallel\pm}=\frac{\omega v_\parallel}{v^2}+\left[\frac{v_\perp^2}{v^2}\left(k_\pm^2-\frac{\omega^2}{v^2}\right)\right]^{1/2}.
$$

前段の式(2)–(12)と後段の式(18)–(20)の転記は `work/transcription/` にある。二次資料と照合する場合も、この表の印刷形を基準にする。

### 4.2 項ごとの対応

項記号は本パッケージで付けたものである。「正孔」の列は、正孔の因子 πδ(ω⊥) + i/ω⊥ を含むかどうかを示す（原著 p. 1615左: 正孔に関係する項はすべてこの因子を含む）。

実装は `src/xpsloss/xps_fixed_depth.py` にある。`_term_integrands` が全8項を評価する。独立経路 I3 との照合は `k_parallel` の分散0・1で行った。分散0では両指定が一致する。

| 項記号 | 式と部分 | 誘電関数の因子 | 正孔 | 1/x | 原著の説明 | 実装 | 検証 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 原著式(1)（[EQUATIONS.md](EQUATIONS.md) の E2） | 式(1) | Im[1/ε]（印刷のまま。Q2） | − | − | 無限媒質 | `model.bulk_diimfp`（M0。Im[−1/ε] で実装） | M0の解析解・和則テスト |
| 13a | 式(13)第1項 | 1/ε − 1 | − | なし | a → ∞ で残る唯一の項（p. 1615左） | `term_13a` | I1、I2 |
| 13b-A/B/C | 式(13)第2項の波括弧の3項 | 1/ε − 1 | − | あり | 媒質内の項（p. 1615左） | `term_13b`（A/B/C の和） | I1、I3 |
| 14-1/2 | 式(14)の波括弧の2項 | 1/ε − 1 | あり | あり | 媒質内、正孔の項 | `term_14`（2項の和） | I3、I5 |
| 15-T1 | 式(15)角括弧の第1項 | (ε−1)(ε+2)/(ε(ε+1)) | − | あり | 媒質内 | `term_15_T1` | I3 |
| 15-T2 | 式(15)角括弧の第2項 | 同上 | あり | あり | 媒質内、正孔の項 | `term_15_T2` | I3、I5 |
| 15-T3 | 式(15)最後の項 | +(ε−1)/(ε+1)（Q4） | − | あり | 媒質内 | `term_15_T3` | I3 |
| 16-e | 式(16)波括弧の分数 | (1−ε)/(1+ε) | − | あり | 真空中の損失（p. 1615左） | `term_16_e` | I3 |
| 16-h | 式(16)波括弧の最後の項 | 同上 | あり | あり | 真空中、正孔の項 | `term_16_h` | I3、I5 |
| E17 | 式(17) | − | − | − | 13a–16-h の積分範囲 | `_parallel_limits` | I2（θ = 0） |
| E18 | 式(18) 経路長の重み付き平均 | − | − | − | 深さ平均 | 本体は未実装（M3）。図8の検算には `benchmarks/paper_figures/model.py` の `depth_average` を使う | − |
| E19 | 式(19) 1回散乱のスペクトル | − | − | − | スペクトル | M4（未実装） | − |
| E20 | 式(20) 原著のAl誘電関数 | − | − | − | 原著の計算例の入力 | 本体の振動子模型の対象外。図3・5・8の検算には `benchmarks/paper_figures/model.py` の `eq20` を使う | − |

M2-1の実装（`xps_fixed_depth.py`）は、ギャップのない振動子の複素応答 E4b を、運動量 q の配列について直接評価する。これが M1 の `dielectric.negative_inverse_dielectric_function` と一致することは試験で確かめている。誘電因子は S = 1 − 1/ε から計算する。1/ε − 1 = −S、(ε−1)/(ε+1) = S/(2 − S)、(ε−1)(ε+2)/(ε(ε+1)) = 2S − S/(2 − S) である。ギャップのある材料の固定深さ計算（M2-2）では、係数を F = −1/ε で作り、節点の実部を残存帯域の閉形式 D3 で求める。ギャップ端と実軸上の極の扱い（遅延極限の点値、広がりなし、ギャップ下の ε = −1 の寄与がないことを証明できる損失だけを計算）は [M2-2 の設計](M2_2_DESIGN.md) で決めた。

「原著の説明」列の「媒質内」「真空中」は原著の区分（p. 1615左）である。各項をどの電荷源（電子、鏡像、正孔、表面電荷）に対応させるかは、2つの転記がそれぞれ導出を示しているが、原著には書かれていない。実装で項の名前に使う場合は「導出による」と記録する。

## 5. 未確定の点

両転記の不確定事項（Fable U1–U12、Codex U1–U9）を統合した。M2-1（θ = 0、ħω > 0）への影響で分類する。

### 5.1 M2-1の結果に影響する

| 番号 | 内容 | 位置 | 扱い |
| --- | --- | --- | --- |
| Q1 | 式(13)–(16)の ε の運動量引数。原著は式(7)の要請として ε(**k**,ω) = ε(k∥,ω) と書き、その妥当性は以前に議論したとして参考文献5を挙げる（p. 1614右）。一方、θ = 0、a → ∞ で13a項が式(1)と一致するという記述（p. 1615左）は、分散があるときは εを k = (k∥² + Ω0²)^{1/2} で評価した場合にだけ成り立つ。数値計算でどちらを使ったかは原著に書かれていない | p. 1614右、p. 1615左 | 式レビューAの結果: 印刷された閉形式は ε(k∥) のときの式(2)の値（分散ありでも数値経路と 10⁻¹⁴ で一致）。分散があるとき、法線出射で13aが式(1)に一致するのは `k_total` である。無分散なら両指定が一致する。M2-1のAl（分散1）では両者で `total` が最大約5%違う。実装では明示的な選択肢とした（[M2_FIRST_CASE.md](M2_FIRST_CASE.md) §5）。参照ソフトとの比較では、公表された式のままだと、k_parallel と k_total のどちらも不合格だった。試験した条件の参照ソフトの出力は、全項で k_total を使い、逆誘電関数の実部の分子に分散前のエネルギーを使う候補の式と一致した（[QUEELS_DISPERSION_TEST.md](QUEELS_DISPERSION_TEST.md)） |
| Q6 | ħω = 0 での δ項と i/ω⊥ の扱い（主値など） | 式(9)、(14)–(16) | ħω > 0 だけを計算する。δ項は寄与しない。原著 p. 1618 の 0.05 eV は式(19)の積分下限であり、K_eff の下限ではない |
### 5.2 解決済み、または結果に影響しない

| 番号 | 内容 | 扱い |
| --- | --- | --- |
| Q4 | 式(15)最後の係数 +(ε−1)/(ε+1) | 印刷どおり正しい（式レビューA、導出と数値の両方）。式(11)の (1−ε)/(1+ε) に、式(10)の F^V の負号が掛かってこの符号になる |
| Q5 | 式(9)は δ(ω)、i/ω、式(14)–(16)は δ(ω⊥)、i/ω⊥ | 印刷どおり正しい（式レビューA）。πδ(ω) + i/ω = (πδ(ω⊥) + i/ω⊥)/v⊥ で、1/v⊥ は共通の前因子に入っている |
| Q14 | 式(2)の前因子 2/((2π)⁴ x ħ² ω) と式(13)–(16)の前因子の整合 | 整合する（式レビューA）。a0 = ħ²/(me²) を使えば厳密に一致する |
| Q2 | 式(1)は Im[1/ε] と負号なしで印刷される。式(20)（Im ε > 0）では負になり、13a項の Re[i(1/ε−1)] = −Im(1/ε) と逆符号である。manual の式(2)、(3)も同じく負号なしで印刷される | 物理量は非負の Im[−1/ε] とする。M0の実装と13a項の印刷形に合わせる。転記は印刷のまま残す |
| Q3 | 式(17)は平方根の前に + だけが印刷される | θ = 0 では k∥± = (k±² − ω²/v²)^{1/2} で、根の選択は生じない。θ > 0 はM2-1の対象外 |
| Q7 | 式(6)第3行の左辺が ρ^V(**r**,t) と印刷される（ρ^V(**k**,ω) の誤植と見られる） | 式(13)–(16)には影響しない |
| Q8 | 速度の字形が ν に似る | 本文の定義により速度 |
| Q9 | 単位系の名称が書かれていない | a₀ を Å、E⊥ を eV とすれば、前因子の単位は eV⁻¹ Å⁻¹ になる（式レビューA5で確認） |
| Q10 | K_eff の引数の並び（(E0,ħω,x)、(E0,ħω,a,θ)、(E0,ħω,θ,a)）が不統一 | 意味は同じ |
| Q11 | 正孔なし・θ = 0 で一致するとされる「模型 A+」（参考文献4）が本文で定義されない | 参考文献4が入手できれば外部照合に使う |
| Q12 | p. 1618の本文で K_sc の引数が式(19)と違う | 式(19)を正とする。M4で扱う |
| Q13 | 図10の凡例が「eq.19」、図の説明は式(20) | 図の表記だけの問題 |
| Q15 | θ > 0 での方位角。d³k = 2πk∥dk∥dk⊥ と書きながら、閉形式は k∥·v∥ を k∥v∥ に置いた（方位角0の）形である（p. 1614） | 式レビューA4の指摘。M2-1（θ = 0）には影響しない。a → ∞ の13a項は、印刷形と式(17)では θ = 60° で式(1)の 0.70–0.76 倍になり、方位角を積分すると θ によらず式(1)に一致した。原著が斜出射の不一致を k⊥ 範囲の近似のためかもしれないとした点（p. 1615左）は、これでも説明できる（推論）。θ > 0 を実装する前に検討する |

## 6. 用語の注意

原著の p. 1615左にある K_eff = K_intr/x + K_extr は、経路長への依存で分けた分解である。図5（p. 1616右）の electron contribution は正孔電荷0の全量である。hole contribution は電荷1と0の全量の差で、干渉を含む。manual では図5の分け方に intrinsic/extrinsic の名を使う。本パッケージの出力名は [M2_FIRST_CASE.md](M2_FIRST_CASE.md) §2 の `total`、`no_hole`、`hole` とする。
