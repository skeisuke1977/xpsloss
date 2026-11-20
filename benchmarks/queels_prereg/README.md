# 分散の処方を判別する事前登録試験

このベンチマークは合成入力 Alg075、Ala025、Al2osc を使い、逆誘電関数の実部の分子を光学極限のエネルギーで評価する候補と、公表式を比較する。候補はこの試験専用である。振動子、参照計算の全設定、式と判定条件は事前登録の「試す仮説」「新しい計算」「合否の基準（固定）」にある。[事前登録](../../docs/QUEELS_DISPERSION_TEST.md)

ソース配布のルートで、同じ版の xpsloss をインストールするか、`PYTHONPATH` にその `src` を指定する。

```sh
python benchmarks/queels_prereg/evaluate_prereg.py predict
python benchmarks/queels_prereg/evaluate_prereg.py predict --out nominal.txt
python benchmarks/queels_prereg/evaluate_prereg.py predict --check
python benchmarks/queels_prereg/evaluate_prereg.py evaluate PATH_TO_YOUR_EXPORTS
python -m unittest discover -s tests -p test_queels_prereg.py
```

`predict` は深さ 10 Å の全予測を ASCII・CRLF で標準出力に書く。Windows PowerShell 5.1 の `>` はテキストを UTF-16 に変換するため、ファイル保存には `--out` を使う。`--out nominal.txt` は同名のファイルを上書きするため、既存の予測を残すには別名を指定する。`--check` の JSON は、再計算したバイト列と同梱ファイルの一致を `identical`、同梱ファイルの SHA-256 と事前登録値の一致を `matches_registered_sha256` で示す。`sha256` は再計算値、`frozen_sha256` は同梱ファイルのハッシュである。両方の照合が通れば終了0、どちらかが不一致なら1。

`evaluate` は指定ディレクトリ直下に、QUEELS-XPS v2.1 で作成した `K_XPS_{T,e,h}_{入力名}_1000_0_q1_n2000.m` を計9個必要とする。ギャップ0、法線から0°、正孔電荷1、1000 eV での実行を使う。`@` 以降のヘッダのエネルギー・角度・IMFP・刻み・点数・ac/bc と、65深さ列の格子を検査する。先頭3行は深さ・固体内の経路長・深さごとの量、その後の500行の損失は行番号 i（1始まり）×0.2 eV で決まる。独立な損失軸の列はない。最終行は未計算のゼロ行。成分間で深さと経路長を照合し、実行ごとの格子を使う。

各実行で10 Åに最も近い列を選び、5.0–99.8 eVの475点を使う。成分 T、e、h の最大絶対差を同じ列・区間の参照 T の最大値で割る。候補の3指標がすべて 0.003 以下、公表式の少なくとも1指標が 0.003 超、候補400/800パネル間の全成分の最大差が同じ分母の 0.0001 以下で当該実行が合格する。3実行すべての合格を要する。

登録文書は「最も浅い列」を診断に指定した。凍結した評価コードが実際に使った添字1（0始まり）は、浅い方から2番目の列である。配布版も浅部に添字1、深部に添字64を使う。この違いは診断値にだけ関係し、合否は各実行で選んだ10 Å最近傍の列から決まる。

標準出力は JSON。終了コード0は全実行が合格、1は仮説未確認、2はファイル不備・数値計算エラーを示す。`runs` 以下の `criterion` に成分別の指標と `candidate_pass`、`discriminating`、`selfcheck_pass` が入る。`diagnostic_shallow` と `diagnostic_deep` は合否に加算しない。

## ファイルの来歴

| ファイル | 事前登録の SHA-256 | 配布版の SHA-256 |
|---|---|---|
| numerator_rule.py | `a57b11d1d8eb0f5000da123389574bc82e6d2d0481a297b69919f80f28f81ccd` | `848dca6b75db80813a8f1ef71a747ea1fb8e59b411061cccebb951a50ad1c6e8` |
| evaluate_prereg.py | `ba213db1344c46b64d511dee004d98b7f8af4219f326f793cbd0a26196e2d77c` | `54c838d7fedf4a4726c4d0756b1e86b267111a63fc16539824dfedc5237b28d2` |
| predictions_a10.txt | `715f4f455227beacce8880f5697a071a531b3f07fd7d93b4eacc8fd5163f2a3d` | `715f4f455227beacce8880f5697a071a531b3f07fd7d93b4eacc8fd5163f2a3d` |

候補コードから非公開の参照値表、それを使った旧実演、不要になった import を除き、説明文を改めた。`terms` と `combine` は同じ演算を同じ順で行う。`terms` の式(17)の呼び出し先を同梱した補助関数に替えた。この関数の本体は登録直後の M2-1 コミットの `xps_fixed_depth._parallel_limits` と共通で、定数名を `H` にした。登録時は未コミットの作業ツリーを使っており、名目予測15行の再現をバイト単位で確かめた。公表式は現行のインストールまたは `PYTHONPATH` のパッケージから `effective_cross_section_fixed_depth` を `k_total` で呼ぶ。
