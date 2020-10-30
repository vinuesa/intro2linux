# Taller de introducci&oacute;n al biocómputo en sistemas GNU/Linux

¡Bienvenid@s al Taller de introducción al biocómputo en sistemas GNU/Linux!

## Justificación y objetivos del taller

El trabajo en gen&oacute;mica se realiza en servidores UNIX o GNU/Linux de alto rendimiento. Recomiendo por tanto
que te familiarices con este ambiente de c&oacute;mputo al inicio de tu formaci&oacute;n acad&eacute;mica. 

Este taller tiene por objetivos iniciarte en el camino del bioc&oacute;mputo en sistemas GNU/Linux y ayudarte a descubrir un ambiente de c&oacute;mputo mucho m&aacute;s amigable y poderoso que el que posiblemente conoces hasta ahora.

Aprender&aacute;s todo lo necesario para un arranque r&aacute;pido, &uacute;til, entretenido y exitoso de programaci&oacute;n/scripting del Linux Shell, usando ejemplos relevantes para el an&aacute;lisis de secuencias moleculares y datos asociados.

## Ediciones del Taller

1a. Edición: semestre 2021-1

v.2020-10-25

***
 
# Presentaci&oacute;n

## El profesor
Hola, me llamo [Pablo Vinuesa](http://www.ccg.unam.mx/~vinuesa/). Soy investigador titular del 
[Centro de Ciencias Gen&oacute;micas](http://www.ccg.unam.mx) de la 
[Universidad Nacional Aut&oacute;noma de M&eacute;xico - UNAM](http://www.unam.mx/).

Mis [l&iacute;neas de investigaci&oacute;n](http://www.ccg.unam.mx/~vinuesa/research.html) 
integran la gen&oacute;mica y la bioinform&aacute;tica con la biolog&iacute;a y gen&eacute;tica molecular para entender 
la evoluci&oacute;n y emergencia de pat&oacute;genos oportunistas a partir de microbios ambientales.

## Sobre el material did&aacute;ctico
A trav&eacute;s de estas p&aacute;ginas se distribuyen los apuntes, ejercicios y datos que se usar&aacute;n en el Taller.
Es una recopilaci&oacute;n del material que he desarrollado para diversos cursos y talleres que ha impartido principalmente en la [Universidad Nacional Aut&oacute;noma de M&eacute;xico - UNAM](https://www.unam.mx/), pero también en otras universidades y países, como Argentina, Brasil, España y Puerto Rico, enlistados en esta p&aacute;gina: [Talleres y Cursos impartidos por Pablo Vinuesa](https://www.ccg.unam.mx/~vinuesa/cursos.html). 


### Licencia y términos de uso
El material docente del curso [**intro2linux**](https://github.com/vinuesa/intro2linux) lo distribuyo p&uacute;blicamente a trav&eacute;s de este [repositorio GitHub](https://github.com/vinuesa/intro2linux) bajo la [**Licencia No Comercial Creative Commons 4.0**](https://creativecommons.org/licenses/by-nc/4.0/) 

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 

<a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html"><img alt="GNU General Public License v3" style="border-width:0" src="https://www.gnu.org/graphics/gplv3-127x51.png" /></a><br />
El código asociado a este Taller se distribuye bajo la licencia <a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html">GNU General Public License v3</a>

### Clonaci&oacute;n del repositorio
Si tienes instalado [git](https://git-scm.com/) en tu computadora, puedes clonar el repositorio con el comando:

   <code>git clone https://github.com/vinuesa/intro2linux.git</code>

- En [ubuntu](https://www.ubuntu.com/) es muy f&aacute;cil instalar git: 

  <code>sudo apt install git</code>

- También en [mobaXterm](https://mobaxterm.mobatek.net/) es posible instalar git como un <i>plugin</i> [mobaXterm - plugins](https://mobaxterm.mobatek.net/plugins.html)

abre una consola de moba (ver instrucciones de instalación abajo) y teclea

  <code>sudo apt-get install git</code>


#### Actualización del repositorio local
Una vez clonado, puedes actualizar tu copia del repo entrando al directorio donde lo clonaste con el comando anterior y ejecutando el siguiente:

<code>git pull</code>

<!--### ¿Horario y lugar de impartici&oacute;n de las sesiones?
Las clases se imparten de manera remota, v&iacute;a zoom, los miércoles de 9:00 - 11:00-->

<!--<img src="docs/pics/intro2linux_aula_UNLP_2-6Julio2018.jpg" />-->


#### Descarga de archivos individuales o del repositorio como archivo zip

Tambi&eacute;n puedes descargar la distribuci&oacute;n como archivo comprimido "zip", o descargar los archivos individuales que desees
desde las carpetas docs y data.

### Instalación de mobaXterm en máquinas windows para poder establecer sesiones ssh a un servidor
Les comparto unas notas elaboradas por la UATI de la @lcg_unam sobre instalación de [mobaXterm - home edition](https://mobaxterm.mobatek.net/download-home-edition.html) en m&aacute;quinas Windows para poder establecer conexiones remotas a servidores v&iacute;a SSH y tener acceso a una consola Linux corriendo <i>bash</i>

- [instalaci&oacute;n mobaXterm - PDF](https://github.com/vinuesa/intro2linux/tree/master/docs/ConexionSSHdesdeWindows_usando_mobaXterm_UATI_LCG-UNAM.pdf)



# Contenidos del Taller de Introducci&oacute;n al bioc&oacute;mputo en sistemas GNU/Linux
En este tutoral se describen los siguientes conceptos b&aacute;sicos:
- Qu&eacute; es el bioc&oacute;mputo?
- Qu&eacute; es GNU/Linux?
- Qu&eacute; es el Shell?
- C&oacute;mo me conecto a un servidor remoto v&iacute;a ssh?
- Revisaremos comandos b&acute;sicos y c&oacute;mo moverte por el sistema de archivos

Si nunca has trabajado en un sistema Linux anteriormente, recomiendo leer este tutoral antes de inciar las sesiones prácticas del siguiente rubro.

- [Primer contacto con Linux - PDF](https://github.com/vinuesa/intro2linux/tree/master/docs/intro_biocomputo_Linux.pdf)

## Tutoral extenso de comandos GNU/Linux con aplicaciones a biocómputo
- [Tutoral - html](https://vinuesa.github.io/intro2linux/docs/)

<!--
#### Pr&aacute;ctica 2. Descarga de secuencias en formato FASTA de GenBank usando el sistema ENTREZ y parseo de los archivos usando herrramientas de filtrado
- [pr&aacute;ctica2 - html](https://vinuesa.github.io/intro2linux/practica2_parseo_fastas/)
- [pr&aacute;ctica2 - pdf](https://vinuesa.github.io/intro2linux/practica2_parseo_fastas/ejercicio_parseo_fastas_ENTREZ.pdf)
- [pr&aacute;ctica2 - fasta](https://vinuesa.github.io/intro2linux/practica2_parseo_fastas/data/recA_Bradyrhizobium_vinuesa.fa)
-->

