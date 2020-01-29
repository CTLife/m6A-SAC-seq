cat nonDRACH.txt | grep -P "\-(100)%\-\S::"        | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(9\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(8\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(7\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(6\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(5\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(4\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(3\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(2\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(1\d[\d.]*)%\-\S::"  | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-([1-9])%\-\S::"      | wc -l >> nonDRACH.number.txt
cat nonDRACH.txt | grep -P "\-(\d\.\d*)%\-\S::"    | wc -l >> nonDRACH.number.txt


