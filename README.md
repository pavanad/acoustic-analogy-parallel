Implementação numérica paralela da analogia acústica
====================================================

Author: Adilson Pavan	RA: 159018  
Unicamp MO644/MC900 (Parallel Programming) 2014s1  

Implementação numérica paralela da analogia acústica de Ffowcs Williams &amp; Hawkings.  


Entrada de dados
----------------

O programa utiliza 3 arquivos para entrada de dados, porém esses dados definem a superfície do aerofólio NACA0012
e o fluxo de pressão em torno do aerofólio. Considerando que o objetivo é calcular o ruído gerado pela pressão ao 
redor do aerofólio, apenas o número de observadores deve ser alterado, sendo assim a aplicação foi desenvolvida
para receber parâmetros que definem as configurações de simulação.  

NACA0012.dat  
pressure_imaginary.in  
real_pressure.in  

Os parâmetros para execução do programa determinam o tamanho e o desempenho do paralelismo.  


| -t        | número de threads (default=1) |
| -o        | número de observadores (default=360)  |
| -serial   | ativa a versão serial (default=false) |
| -parallel | ativa a versão paralela (default=true) |  
| -output   | ativa a saída no padrão da disciplina MO644/MC900 (default=true) |

exemplo: ./surf -t 4 -o 360 -serial -parallel -output input01  


Resultados
----------

Além do arquivo output.dat gerado no padrão da disciplina os resultados estão no padrão do software Mathematica para facilitar a apresentação dos dados.  

| NACA0012.nb				| Representação do aerofólio NACA0012 no padrão Mathematica |
| PressaoAtTeta_input01.nb	| Ruído gerado pela pressão ao redor do aerofólio |
| PressaoAtPoint_input01.nb	| Os dados de entrada de pressão no padrão Mathematica |
