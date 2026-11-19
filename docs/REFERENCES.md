# 一次資料と用途

確認日: 2026-11-19。PDF本文や元ソフトは本パッケージに同梱していない。

## R1. 公式製品ページ

QUEELS-XPS, Sven Tougaard and Francisco Yubero.

https://www.quases.com/products/queels-xps/

機能の範囲、公式マニュアル、ソフト配布、非商用無料の案内を確認した。公開資料からCLI/API仕様やソース公開は確認できなかった。

## R2. 基礎理論

A. Cohen Simonsen, F. Yubero, S. Tougaard, “Quantitative model of electron energy loss in XPS”, Physical Review B 56, 1612-1619 (1997).

https://doi.org/10.1103/PhysRevB.56.1612

著者公開の本文: https://www.academia.edu/8468514/Quantitative_model_of_electron_energy_loss_in_XPS

バルクの式(1)、半無限媒質と正孔を含む式(13)-(17)、深さ平均・スペクトルへの接続が重要。数式抽出には崩れがあるため、本体の実装前に原PDF画像との照合が必要。

## R3. マニュアル Version 2.1

S. Tougaard, F. Yubero, “QUEELS-XPS: Quantitative Analysis of Electron Energy Loss in XPS Spectra”, User's Guide, Version 2.1。論文庫キー `Tougaard2011` の年はp.2の著作権表示（2001–2011）の末年に基づく。刊行日の明記はない。

https://www.quases.com/data/documents/Ver=202.1=20Manual=20QUEELS-XPS-simulation.pdf

理論の概要、入力・出力、Si例、GUI操作、フォーマット、モデルの制限。Si入力はp.26表1、計算条件はp.27 §8.1.2、角度の基準はp.7図1およびp.20 §4.3。論文庫のPDFは[R9]から取得した。

## R4. QUEELS全体の方法論

S. Tougaard, N. Pauly, F. Yubero, “QUEELS: Software to calculate the energy loss processes in TEELS, REELS, XPS and AES including effects of the core hole”, Surface and Interface Analysis 54, 820-833 (2022).

https://doi.org/10.1002/sia.7095

https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/sia.7095

バルクDIIMFPは式(12)。各幾何の違い、直線軌道の仮定、ELFの役割、参照例を確認。

## R5. 一次励起を取り出す手順とCuO

N. Pauly, F. Yubero, S. Tougaard, “Determination of the Primary Excitation Spectra in XPS and AES”, Nanomaterials 13, 339 (2023).

https://doi.org/10.3390/nano13020339

https://www.mdpi.com/2079-4991/13/2/339

ELFの式(3a),(3b)、一次励起と損失の分担、CuOの電荷移動多重項計算との比較を確認。

## R6. Alの角度依存

F. Yubero, S. Tougaard, “Quantification of plasmon excitations in core-level photoemission”, Physical Review B 71, 045414 (2005).

https://doi.org/10.1103/PhysRevB.71.045414

論文庫キー `Yubero2005`。Al 2sの角度依存を、表面・バルク損失の定量的な参照にする。p.045414-3に記載されたAl振動子は[R8]のAl入力と異なる。

## R7. NiOと非局所遮蔽

N. Pauly, F. Yubero, F. J. García-García, S. Tougaard, “Quantitative analysis of Ni 2p photoemission in NiO and Ni diluted in a SiO2 matrix”, Surface Science 644, 46-52 (2016).

https://doi.org/10.1016/j.susc.2015.09.012

著者所属機関の概要: https://portal.findresearcher.sdu.dk/en/publications/quantitative-analysis-of-ni-2p-photoemission-in-nio-and-ni-dilute/

誘電損失と局所・非局所終状態を組み合わせた解析の先例。CoOに同じ結論が成立するという根拠ではない。

## R8. 材料のELFデータ集

N. Pauly, F. Yubero, S. Tougaard, “ELF dielectric functions for various materials”, Zenodo, Version 1.0 (2022-02-09), CC-BY-4.0.

https://doi.org/10.5281/zenodo.6024064

https://zenodo.org/records/6024064

`benchmarks/sources/zenodo6024064_v1.0/` に `ELF_Si.txt`、`ELF_Au.txt`、`ELF.txt` の原本を収録（2026-11-19取得）。`ELF.txt` のAlと単独ファイルのSi・Auを固定入力として採用した。`ELF.txt` にCoOの収録はない。帰属表示と原典内の行番号は[固定入力と出典](BENCHMARK_INPUTS.md)に示す。

## R9. 2026年の配布資料

S. Tougaard, F. Yubero, “QUEELS-XPS: Software to calculate the energy loss processes in XPS and AES including the effects of the core hole”, Zenodo record v2, published 2026-10-16.

https://doi.org/10.5281/zenodo.23414559

https://zenodo.org/records/23414559

インストールZIP、Version 2.1のmanual、注意書き、関連論文を掲載。登録のv2とソフト本体の新版を同一視しない。

論文庫の `Tougaard2011`、`Yubero2005`、`Yubero2012`、`Tougaard2004`、`Pauly2010` のPDFはこの登録から取得した。本リポジトリはPDFを再配布しない。

## R10. 初期ソフト紹介

S. Tougaard, F. Yubero, “Software package to calculate the effects of the core hole and surface excitations on XPS and AES”, Surface and Interface Analysis 44, 1114-1118 (2012).

https://doi.org/10.1002/sia.4855

干渉を含む一段階のモデルと、ソフトの適用範囲を説明。

## R11. 内殻光電子の損失構造と角度

F. Yubero, S. Tougaard, “Dielectric description of the angular dependence of the loss structure in core level photoemission”, Journal of Electron Spectroscopy and Related Phenomena 185, 552–558 (2012)。

論文庫キー `Yubero2012`。https://doi.org/10.1016/j.elspec.2012.11.011

## R12. QUEELSによる表面効果の計算

S. Tougaard, F. Yubero, “QUEELS software package for calculation of surface effects in electron spectra”, Surface and Interface Analysis 36, 824–827 (2004)。

論文庫キー `Tougaard2004`。https://doi.org/10.1002/sia.1774

## R13. XPSの表面・内殻正孔効果

N. Pauly, S. Tougaard, “Surface and core hole effects in X-ray photoelectron spectroscopy”, Surface Science 604, 1193–1196 (2010)。

論文庫キー `Pauly2010`。https://doi.org/10.1016/j.susc.2010.04.001

## R14. 逆誘電関数のクラマース・クローニッヒ関係

V. B. Bobrov, S. A. Trigger, G. J. F. van Heijst, P. P. J. M. Schram, “Kramers–Kronig relations for the dielectric function and the static conductivity of Coulomb systems”, EPL 90, 10003 (2010).

https://doi.org/10.1209/0295-5075/90/10003

著者公開版: https://arxiv.org/abs/1003.4724

[EQUATIONS.md](EQUATIONS.md) E4と [DIELECTRIC_DESIGN.md](DIELECTRIC_DESIGN.md) が、逆誘電関数の分散関係と周波数対称性の根拠として著者公開版 p. 2 の式(1),(2),(4),(5),(7)を引用する。論文庫には未収録で、式番号は実装担当の参照に基づき、独立レビューでは照合していない。E4aの符号・係数・高周波極限そのものは、独立レビューで導出し直して確認した（[VALIDATION.md](VALIDATION.md)）。
