# M1-B 振動子 ELF からの複素誘電応答

損失エネルギー \(T=\hbar\omega\) と積分変数 \(u\) は eV、波数 \(q\) は Å\(^{-1}\) とする。誘電関数 \(\epsilon\)、\(F=-1/\epsilon\)、正損失での ELF \(L=\operatorname{Im}F\) は無次元。固定した \(q\) で振動子中心を \(\Omega_i=E_i+\alpha_i\hbar^2 q^2/(2m_e)\)、ギャップを掛ける前の ELF を \(L_0\) とする。

## 符号と幅

QUEELS-XPS v2.1 マニュアル p. 10 の式 (2) と (3) は、どちらも \(\operatorname{Im}\{1/\epsilon\}\) と印刷されている。式 (3) の振動子和は \(T>0\)、非負強度で非負になる。Tougaard ら (2022) PDF p. 4 の式 (9) と p. 5 の式 (12) は、同じ正の ELF を \(\operatorname{Im}\{-1/\epsilon\}\) と記す。一方、同じ p. 5 の式 (11) は左辺を \(\operatorname{Im}\{1/\epsilon\}\) と印刷している。以下では \(e^{-i\omega t}\) を採用する。この規約では受動応答の \(\operatorname{Im}\epsilon\) と \(L\) は非負で、\(\operatorname{Im}\{1/\epsilon\}=-L\) となる。式 (11) の右辺を保つには式 (9) と同じ \(L\) を左辺に入れる必要がある。マニュアルの式 (2)、(3) にもこの符号を対応させる。マニュアルの印刷がフーリエ規約の違いに由来するかは、該当頁からは確定できない。[マニュアル v2.1][manual]、[Tougaard ら (2022)][tougaard]

マニュアル p. 26 の Table 1 は \(\gamma_i\) を “half width of peaks” と呼ぶ。p. 10 の式 (3) は \(\gamma_i\hbar\omega\) をそのまま二乗項に使う。`width_eV` はこの分母に入る減衰パラメータ \(\Gamma_i\) であり、表の値を式 (3) のパラメータとして転記する場合は係数 2 を掛けない。狭い単一ピークでは \(|T-\Omega_i|,\Gamma_i\ll\Omega_i\) として

\[
L_{0,i}(T)\simeq\frac{A_i\Gamma_i}
{\Omega_i[4(T-\Omega_i)^2+\Gamma_i^2]}.
\tag{D1}
\]

この近似で半値半幅は \(\Gamma_i/2\)、半値全幅は \(\Gamma_i\) になる。広いピークやギャップに近いピークでは、式 (9) の線形状自体から幅を求める。GUI 内部の幅変換は未照合。[マニュアル v2.1][manual]、[Tougaard ら (2022)][tougaard]

## 振動子入力

ギャップがゼロなら [E4b](EQUATIONS.md#e4-振動子-elf-からの複素誘電応答) を使う。各極は \(z=-i\Gamma_i/2\pm\sqrt{\Omega_i^2-\Gamma_i^2/4}\) で、\(\Omega_i,\Gamma_i>0\) なら下半平面にある。\(F(\infty)=-1\) と \(\epsilon(\infty)=1\) を満たす。\(T=0\) の静的逆誘電関数は \(1-\sum_i A_i/\Omega_i^2\) である。

ギャップが正なら、実部も切断後の ELF に合わせて E4a のクラマース・クローニッヒ（KK）変換で求める。Bobrov らの式 (7) は \(1/\epsilon\) の高周波定数 1 と Cauchy 主値を含む。両辺を負にし、式 (4) の負周波数対称性を \(1/\epsilon\) にも適用すると、正のエネルギー上の核 \(2u/(u^2-T^2)\) と定数 \(-1\) が得られる。式 (1)、(2) の遅延応答と式 (5) の正の \(\operatorname{Im}\epsilon\) も規約と整合する。ギャップの位置に応じて、除去帯域を解析応答から引くか残存帯域を積分する数値経路を標準にする。全帯域を振動子の形で表すため、出力格子の上限で積分は打ち切らない。[Bobrov ら (2010)][bobrov]、[著者公開版の PDF p. 2][bobrov-pdf]

### 独立な部分分数評価

一つの振動子を取り、添字 \(i\) を省く。\(g=E_g>0\)、\(T\ne g\) とし、
\(D(u)=(\Omega^2-u^2)^2+\Gamma^2u^2\) とおく。式 E4a の有理関数は

\[
\frac{u^2}{(u^2-T^2)D(u)}
=\frac{1}{D(T)}\left[
\frac{T^2}{u^2-T^2}
+\frac{\Omega^4-T^2u^2}{D(u)}\right].
\tag{D2}
\]

\(p,b=\Gamma/2\pm\sqrt{\Gamma^2/4-\Omega^2}\) とすると、\(p+b=\Gamma\)、\(pb=\Omega^2\)、\(D(u)=(u^2+p^2)(u^2+b^2)\)。主値を含む原始積分と、因数分解した残りを使うと、振動子が \(\operatorname{Re}F+1\) に加える量は

\[
\begin{aligned}
P_g(T)&=\operatorname{PV}\int_g^\infty\frac{du}{u^2-T^2}
 =\frac{1}{2T}\log\left|\frac{g+T}{g-T}\right|,\qquad
 P_g(0)=\frac{1}{g},\\
J_g(p)&=\int_g^\infty\frac{du}{u^2+p^2}
 =\frac{\arctan(p/g)}{p},\qquad
Q_g(p,b)=\frac{J_g(p)-J_g(b)}{b^2-p^2},\\
\Delta\operatorname{Re}F
&=\frac{2A\Gamma}{\pi D(T)}\operatorname{Re}\left[
T^2\{P_g(T)-J_g(b)\}+(\Omega^4+T^2p^2)Q_g(p,b)
\right].
\end{aligned}
\tag{D3}
\]

複素 \(\arctan\) は \(\operatorname{Re}p,\operatorname{Re}b>0\) での主値を使う。臨界減衰 \(p=b=\Omega\) では \(Q_g=[J_g(p)-g/(g^2+p^2)]/(2p^2)\)。二根が近い場合は \(\arctan(p/g)-\arctan(b/g)=\arctan[g(p-b)/(g^2+pb)]\) で差を評価する。この式は振動子ごとの解析的な参照経路を与える。狭い振動子がギャップに覆われ、巨大な項の差が残る場合は丸め誤差の指標を確認し、閉形式で解けなければ数値経路を使う。試験では低減衰、臨界減衰の前後、過減衰を含めて、残存帯域を直接 Cauchy 主値積分した参照値と比較する。許容差は \(\text{atol}=10^{-10}\)、\(\text{rtol}=10^{-9}\)。

## 静的な規格化と第一モーメント

光学極限 \(q=0\) で次を診断値として返す。

\[
\begin{aligned}
S_0&=\frac{2}{\pi}\int_{E_g}^\infty \frac{L(0,u)}{u}\,du,
 &R_0&=S_0-(1-r_{\rm ref}),\\
M_1&=\int_{E_g}^\infty uL(0,u)\,du,
 &(S_0,M_1)_{E_g=0}
 &=\left(\sum_i\frac{A_i}{E_i^2},\ \frac{\pi}{2}\sum_i A_i\right).
\end{aligned}
\tag{D4}
\]

\(S_0\) と残差 \(R_0\) は無次元、\(M_1\) は eV²。参照 \(r_{\rm ref}\) は利用者が与える \(\operatorname{Re}\{1/\epsilon(0)\}\) であり、条件が合えば絶縁体の \(1/n^2\) を使える。ギャップゼロの \(S_0\) は E4b の静的極限、\(M_1\) は虚軸上の高エネルギー係数 \(\lim_{\eta\to\infty}\eta^2[F(i\eta)+1]=(2/\pi)M_1\) から得る。ギャップが正なら \(x=E_g/u\) で残存帯域を直接積分する。実材料の電子密度との f 和則の比較には独立した参照が要る。診断結果に基づく入力強度の再規格化は行わない。[Tougaard ら (2022) PDF p. 5、式 (11)][tougaard]

Zenodo 6024064 v1.0 の `ELF.txt` の `Al [4]` は \(E_g=0\)、\(n=999\)、単一振動子 \((E_i,A_i,\Gamma_i,\alpha_i)=(15.6,243.3598,1.5,1)\) を収録する。この行から \(S_0=0.999999178172\)、\(1-1/999^2=0.999998997997\)、残差 \(+1.802\times10^{-7}\) を得た。\(q=0,T=0\) の \(F\) は \(-8.218\times10^{-7}\)、\(\epsilon\) は約 \(1.2168\times10^6\) となる。\(T=0.1,1,15.6\) eV では有限の応答と正の \(\operatorname{Im}\epsilon\) を確認した。厳密な Drude 条件 \(A=E_i^2\) では \(F(0)=0\)。この点で `negative_inverse_dielectric_function` はゼロを返し、`dielectric_response` は極を報告する。正の \(T\) では \(\epsilon=1-A/(T^2+i\Gamma T)\) を評価できる。[Zenodo 6024064 v1.0][zenodo]

### ギャップ端と極

硬いギャップでは、\(L_0(q,E_g)>0\) のとき \(\operatorname{Re}F\) が端の両側で対数発散する。発散係数は \(L_0(q,E_g)/\pi\) で、端そのものの有限値は定義できない。\(F(q,0)<0\) なら、ギャップ下で \(F\) は単調増加して一度ゼロを横切る。この点は \(\epsilon\) の無減衰の実軸極に当たる。これは硬いしきい値を持つモデルの帰結であり、この極を実材料の励起と同定するには別の根拠が要る。実周波数の点評価では極の位置でエラーを返す。極を含む積分には遅延極限と分布としての虚部の扱いが必要になる。\(F(q,0)>0\) なら上半平面に \(F=0\) があるため、`dielectric_response` は因果的な \(\epsilon\) としての評価を拒否する。\(F\) 単体はこの条件でも評価できる。[E4a、E4d](EQUATIONS.md#e4-振動子-elf-からの複素誘電応答)、[Bobrov らの §I][bobrov-pdf]

後段で分母を含む組合せを積分する際は、実軸上の極やギャップ端を通常の格子点として処理せず、原著の遅延規約に従って極限や分布項を確認する必要がある。有限の立ち上がりを採用する場合は ELF のモデルを変更し、その分散実部も再計算する。1997 年の表面・正孔の各項が \(\epsilon\)、\(1/\epsilon\)、\(1/(\epsilon+1)\) のどれを必要とするかは未確定であり、原著の該当式と積分経路の照合を待つ。[Simonsen ら (1997)][simonsen]

## 有限帯域の表形式 ELF

`tabulated_elf_kk` は一つの波数で得た非負の ELF を区分線形で補間し、主値を含む実部を参照計算する。測定区間の両端を \(a,b\) とすると、帯域内の寄与は

\[
\Delta\operatorname{Re} F_{[a,b]}(T)=
\frac{1}{\pi}\operatorname{PV}\int_a^b L_{\rm lin}(u)
\left(\frac{1}{u-T}+\frac{1}{u+T}\right)du.
\tag{D5}
\]

区間内または端の近くでは、Cauchy 核の分子から評価点または近い端の ELF を引く。定数項を対数で計算し、正則化した補間関数を数値積分する。低エネルギー側は `zero` で最初の点までゼロ、または `linear_to_origin` で原点へ直線接続する。\(u=0\) の標本値はゼロとする。先頭が既知のギャップ端で立ち上がりの右極限を与える場合に `zero` を使う。ギャップ端より上まで欠測していれば、最初の点までゼロとする設定は下端の感度評価に相当する。線形接続は低エネルギーで \(L\) が線形に消えるとみなせる場合の外挿である。

高エネルギー側の `hard_cut` は \(u>b\) をゼロとする。`power_law_3` は最後の標本 \(L_b\) に連続に接続する \(L_b(b/u)^3\) を使う。係数 \(L_b b^3\) と末尾二点間の対数傾きを返す。指数は振動子の漸近形から \(-3\) に固定し、末尾の傾きはこの仮定を点検する診断値とする。\(r=T/b\) とした尾部の寄与は

\[
\Delta\operatorname{Re} F_{(b,\infty)}(T)
=\frac{2L_b}{\pi r^2}
\left[\frac{1}{2r}\log\left|\frac{1+r}{1-r}\right|-1\right],
\qquad
\Delta\operatorname{Re} F_{(b,\infty)}(0)=\frac{2L_b}{3\pi}.
\tag{D6}
\]

連続接続では \(T=b\) の両側の対数項を合わせて有限な極限を取る。切断端に跳びがあれば、端そのものは対数発散として報告する。出力には元の標本区間、上下の外挿、補間法、尾部係数と傾き、数値積分設定を含める。数値積分の誤差推定は、欠測域、測定誤差、標本化に起因する不確かさを含まない。

### 上限と刻みに対する感度

合成例は \(q=0,E_g=0,E_i=12\) eV、\(A_i=50\) eV²、\(\Gamma_i=3\) eV。標本は 0 から上限まで等間隔に取り、低エネルギー側は `zero` とした。評価点 \(T=0,3,8,10.25,12,14.75,19\) eV で E4b の解析値と比較した。表は \(\operatorname{Re}F\) の無次元絶対誤差。

| 上限 (eV) | 刻み (eV) | 打ち切り \(T=0\) | 打ち切り最大 | \(u^{-3}\) 尾部 \(T=0\) | \(u^{-3}\) 尾部最大 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 24 | 0.5 | 3.10031e-3 | 8.57382e-3 | 8.82551e-4 | 6.69891e-3 |
| 24 | 0.125 | 3.18861e-3 | 5.70567e-3 | 7.94246e-4 | 1.21518e-3 |
| 24 | 0.03125 | 3.19406e-3 | 5.69463e-3 | 7.88795e-4 | 1.14976e-3 |
| 48 | 0.5 | 2.14505e-4 | 5.48606e-3 | 1.11526e-4 | 5.82489e-3 |
| 48 | 0.125 | 3.04239e-4 | 6.07120e-4 | 2.17908e-5 | 3.39370e-4 |
| 48 | 0.03125 | 3.09779e-4 | 3.44241e-4 | 1.62513e-5 | 3.55238e-5 |
| 96 | 0.5 | 5.90234e-5 | 5.77159e-3 | 9.61152e-5 | 5.80903e-3 |
| 96 | 0.125 | 3.07328e-5 | 3.24968e-4 | 6.35898e-6 | 3.23495e-4 |
| 96 | 0.03125 | 3.62737e-5 | 5.50424e-5 | 8.18084e-7 | 1.96476e-5 |

刻み 0.5 eV では、上限を 48 から 96 eV に伸ばしても最大誤差は約 \(6\times10^{-3}\) に残る。上限 96 eV、刻み 0.125 eV では尾部を付けると静的点の誤差は約 5 分の 1 になるが、評価点全体の最大誤差は同程度である。同じ上限で刻みを 0.03125 eV にすると、尾部を付けた場合の最大誤差は \(1.96\times10^{-5}\) になる。ピーク付近の標本化と測定帯域外の寄与をそれぞれ評価する必要がある。

```python
import numpy as np
from xpsloss.dielectric import tabulated_elf_kk
from xpsloss.model import Material, Oscillator, loss_function

material = Material("synthetic", 0, (Oscillator(12, 50, 3),))
energies = np.array([0, 3, 8, 10.25, 12, 14.75, 19])
reference = (-1 + 50 / (12**2 - energies**2 - 3j * energies)).real
for upper in (24, 48, 96):
    for step in (0.5, 0.125, 0.03125):
        grid = np.linspace(0, upper, round(upper / step) + 1)
        data = loss_function(material, grid)
        for tail in ("hard_cut", "power_law_3"):
            result = tabulated_elf_kk(grid, data, energies, low_energy="zero", high_energy=tail)
            error = abs(result.real_negative_inverse_epsilon - reference)
            print(upper, step, tail, error[0], max(error))
```

### 用途別の選択

| 入力と用途 | 推奨する方法 | 制約 |
| --- | --- | --- |
| ギャップなしの振動子 | E4b の解析式 | 振動子の高エネルギー形と有限 \(q\) の分散は入力モデルに依存 |
| 硬いギャップ付きの振動子 | E4a の数値主値積分を標準、D2–D3 を独立照合に使う | 端や \(\epsilon\) の極を含む点評価に制約。閉形式は強い桁落ちに注意 |
| 有限帯域の測定 ELF | D5 を診断に使い、打ち切りと根拠のある外挿の感度を比較 | 帯域外の強度は標本だけから決まらない |
| 末尾が多数のピークを越えて漸近域にある測定 ELF | 最後の点に接続した \(u^{-3}\) 尾部を試し、末尾の傾きと静的規格化を確認 | 末尾付近に未収録の励起があれば指数だけで補えない |
| 帯域外の独立なデータか物性モデルがある場合 | 測定区間との接続を確認して物性モデルの尾部を使う | 設定の根拠と不確かさを記録する。表形式 API には振動子尾部の接続を追加していない |

## 式、実装、検証の対応

| 式 | 出典と仮定 | 単位 | 実装または確認テスト |
| --- | --- | --- | --- |
| E1、E4b | Tougaard ら (2022) PDF p. 4 式 (9)、(10)。E4b はギャップなしで E1 の虚部を再現する極表示 | \(A_i\) は eV²、\(\Omega_i,\Gamma_i,T\) は eV、\(F,L\) は無次元 | `loss_function`、`negative_inverse_dielectric_function`。単一 Lorentz 振動子、有限 \(q\)、Drude 極限 |
| E2 の ELF の符号 | 2011 マニュアル p. 10 式 (2)、(3) と Tougaard ら (2022) PDF pp. 4–5 式 (9)、(12) を比較 | \(K\) は Å\(^{-1}\) eV\(^{-1}\) | `bulk_diimfp` の解析基準と非負性 |
| D1 | 2011 マニュアル p. 10 式 (3)、p. 26 Table 1。狭い無切断ピークへの展開 | \(L\) は無次元、幅は eV | `width_eV` を係数として使う約束。Lorentz 解析形と ELF 一致テストで係数を確認 |
| E4a、E4c | Bobrov ら著者公開版 PDF p. 2 式 (4)、(7)。周波数の偶奇性と固定 \(q\) の遅延応答。E4c は分子定数の引き算 | \(F,L\) は無次元、\(u,T\) は eV | `negative_inverse_dielectric_function`。直接の Cauchy 主値積分と対数発散係数 |
| D2、D3 | E4a と Tougaard ら (2022) PDF p. 4 式 (9)、(10) から部分分数で導出。ギャップ正、評価点は端と異なる | \(D\) は eV⁴、\(P_g,J_g\) は eV\(^{-1}\)、\(Q_g\) は eV\(^{-3}\) | `negative_inverse_dielectric_closed_form`。減衰 5 条件とギャップ 2 条件の独立主値積分 |
| E4d と極 | \(F=-1/\epsilon\) の定義と E4a の複素延長。Bobrov ら著者公開版 PDF p. 2 式 (1)、(2)、(5)、(7) を参照 | \(F,\epsilon\) は無次元 | `dielectric_response`。真空、受動性、静的な Drude 極、上半平面の極 |
| D4 | Tougaard ら (2022) PDF p. 5 式 (11) の符号を E4a の静的極限から確定。第一モーメントは E4a と E4b の高エネルギー展開 | \(S_0,R_0\) は無次元、\(M_1\) は eV² | `optical_sum_rule_lhs`、`optical_sum_rule_residual`、`optical_elf_first_moment`。ギャップゼロの閉形式とギャップありの直接積分 |
| D5、D6 | Bobrov ら著者公開版 PDF p. 2 式 (7) と E4a の正エネルギー核を、線形補間と固定指数の尾へ適用 | \(a,b,T\) は eV、尾部係数は eV³、\(\Delta\operatorname{Re}F\) は無次元 | `tabulated_elf_kk`。区分線形の解析積分、連続接続の極限、上限と刻みの感度表 |

[manual]: https://zenodo.org/records/23414559
[tougaard]: https://doi.org/10.1002/sia.7095
[bobrov]: https://doi.org/10.1209/0295-5075/90/10003
[bobrov-pdf]: https://arxiv.org/pdf/1003.4724#page=2
[zenodo]: https://doi.org/10.5281/zenodo.6024064
[simonsen]: https://doi.org/10.1103/PhysRevB.56.1612
