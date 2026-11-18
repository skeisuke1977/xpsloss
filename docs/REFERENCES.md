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

https://www.quases.com/data/documents/Ver=202.1=20Manual=20QUEELS-XPS-simulation.pdf

理論の概要、入力・出力、Si例、GUI操作、フォーマット、モデルの制限。特にpp.10-14、20-22、31を参照。

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

Al 2sの角度依存を、表面・バルク損失の定量的な参照にする。

## R7. NiOと非局所遮蔽

N. Pauly, F. Yubero, F. J. García-García, S. Tougaard, “Quantitative analysis of Ni 2p photoemission in NiO and Ni diluted in a SiO2 matrix”, Surface Science 644, 46-52 (2016).

https://doi.org/10.1016/j.susc.2015.09.012

著者所属機関の概要: https://portal.findresearcher.sdu.dk/en/publications/quantitative-analysis-of-ni-2p-photoemission-in-nio-and-ni-dilute/

誘電損失と局所・非局所終状態を組み合わせた解析の先例。CoOに同じ結論が成立するという根拠ではない。

## R8. 材料のELFデータ集

N. Pauly, F. Yubero, S. Tougaard, “ELF dielectric functions for various materials”, Version 1.0 (2022).

https://doi.org/10.5281/zenodo.6024064

https://zenodo.org/records/6024064

Si/Au等の入力と出典。ELF.txtには多材料が含まれるが、今回の確認ではCoOは未収録。データの採用時はソフト本体と区別して利用条件を確認する。

## R9. 2026年の配布資料

S. Tougaard, F. Yubero, “QUEELS-XPS: Software to calculate the energy loss processes in XPS and AES including the effects of the core hole”, Zenodo record v2, published 2026-10-16.

https://doi.org/10.5281/zenodo.23414559

https://zenodo.org/records/23414559

インストールZIP、Version 2.1のmanual、注意書き、関連論文を掲載。登録のv2とソフト本体の新版を同一視しない。

## R10. 初期ソフト紹介

S. Tougaard, F. Yubero, “Software package to calculate the effects of the core hole and surface excitations on XPS and AES”, Surface and Interface Analysis 44, 1114-1118 (2012).

https://doi.org/10.1002/sia.4855

干渉を含む一段階のモデルと、ソフトの適用範囲を説明。
