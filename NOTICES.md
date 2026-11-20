# ライセンスと帰属表示

Keisuke Sato は xpsloss のコードを MIT、自作文書を CC BY 4.0 で提供する。いずれの許諾も下記の第三者データを対象外とする。第三者データには原著者による CC BY 4.0 が適用される。[MITライセンス](LICENSE)・[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

| 区分 | 対象 | ライセンス |
| --- | --- | --- |
| コード・実行例・ビルド設定 | `src/**`、`tests/**`、`benchmarks/**/*.py`、`examples/**`、`pyproject.toml`、`MANIFEST.in`、文書中の自作コード例 | MIT |
| 自作文書と自作の計算結果 | ルートの Markdown ファイル、`docs/**`、`benchmarks/**/*.md`、`benchmarks/queels_prereg/predictions_a10.txt`、`benchmarks/paper_figures/comparison.csv` の計算値、`readings.csv` の読取り値以外の列（図・パネル・条件・読取り誤差の見積り）、`CITATION.cff` の自作部分 | CC BY 4.0。コード例は上記の MIT に従う |
| 原著の図から読み取った値 | `benchmarks/paper_figures/readings.csv` の読取り値と、`comparison.csv`・`RESULTS.md` に再掲した読取り値 | 下の「文書中の引用」のとおり、出典を示した引用として扱い、本プロジェクトの CC BY 4.0 の対象外とする |
| 第三者データ | `benchmarks/sources/**`、`benchmarks/inputs/*elfdb*` | Pauly、Yubero、Tougaard による CC BY 4.0 |

## 文書中の引用

文書には A. Cohen Simonsen, F. Yubero, S. Tougaard, Phys. Rev. B **56**, 1612 (1997) の式(1)、(13)–(17)を検証用に組み直した転記と、図3・5・8から読み取った少数の値と、Table I の線幅5組の転記（`benchmarks/paper_figures/generate.py` と関連文書）、QUEELS-XPS v2.1 User's Guide の短い語句と頁・図表への参照がある。図の画像、切り抜き、重ね描きは同梱しない。引用部分の権利は各権利者に帰属し、出典を示して引用した第三者の素材は本プロジェクトの CC BY 4.0 の対象外である。[式の対応表 §4.1](docs/M1_EQUATION_MAP.md)・[幅の規約と出典](docs/DIELECTRIC_DESIGN.md)

## 同梱データ

N. Pauly, F. Yubero, S. Tougaard, *ELF dielectric functions for various materials*, Version 1.0, Zenodo, 2022-02-09. [データ集](https://doi.org/10.5281/zenodo.6024064)・[原著者のライセンス](https://creativecommons.org/licenses/by/4.0/)

| 場所 | 内容と変更の有無 |
| --- | --- |
| `benchmarks/sources/**` | ELFデータ集の原本3ファイル。原本のバイト列を保存した |
| `benchmarks/inputs/*elfdb*` | データ集のSi・Au・Alから材料パラメータを抜粋してJSONに転記した。来歴には原本の行番号とSHA-256を記録した。強度の再規格化は行っていない。`*.bulk.json` の電子エネルギーと損失格子は本プロジェクトが動作確認用に選んだ |

`examples/` の入力はプロジェクトが作成した合成例である。データ集の固定入力については、原本との照合と、実測に照らした物性値の検証を区別する。実測との比較は未実施である。[固定入力の出典と検算](docs/BENCHMARK_INPUTS.md)

## 参照ソフトとの関係

QUEELS-XPS および QUASES は各権利者のソフトウェアの名称であり、ここでは参照対象を特定するために用いる。

xpsloss は公表された式から独立に実装した。QUEELS-XPSの公式版・承認版・互換保証版ではない。元製品の実行ファイル、DLL、ソースコード、図版、論文PDFは同梱しない。外部データを追加する際は、その資料の利用条件を個別に記録する。

分散のあるAlでは、QUEELS-XPS v2.1 の出力との比較で両方の運動量指定が事前の基準を満たさなかった。続く合成入力3件では、試験した条件の出力が実部の分子に分散前のエネルギーを使う候補式と所定の許容幅で一致した。本パッケージは Yubero ら（1996）の公表式(32)を使う。比較の条件と集約した差を公表し、参照ソフトの出力ファイルと数値表は非公開とする。[当初の比較](docs/M2_FIRST_CASE.md)・[追試の条件と結果](docs/QUEELS_DISPERSION_TEST.md)

## English summary

Keisuke Sato licenses the code (including examples in the documentation) under MIT and the project's own documentation under CC BY 4.0. Quoted equations, terms from manuals, and other cited third-party material retain their owners' terms and are excluded from the project's CC BY 4.0 grant. Both project grants exclude the original data in `benchmarks/sources/**` and their derivatives in `benchmarks/inputs/*elfdb*`, which retain the CC BY 4.0 license of N. Pauly, F. Yubero and S. Tougaard, *ELF dielectric functions for various materials*, v1.0 (2022). xpsloss is independent of QUEELS-XPS, QUASES and their authors and does not distribute reference-software output. [ELF data set](https://doi.org/10.5281/zenodo.6024064) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
