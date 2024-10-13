# How to Git
Git ist ein leistungsstarkes Versionskontrollsystem, mit dem man Änderungen
an Dateien verfolgen, mit anderen zusammenarbeiten und verschiedene Versionen
von Quellcodes verwalten kann.
GitHub ist eine Plattform, auf der Git-Repositories online gehostet werden 
und die eine einfache Zusammenarbeit ermöglicht. Dieses Tutorial geht auf die
Einrichtung von Git, die Zusammenarbeit an Projekten und die Nutzung von
Github ein

## Installieren von Git
Für Debian passierte Distributionen kann man einfach über
```
sudo apt install git-all
```
installieren. Für weitere Betriebssysteme kann man dies auf der [Website](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) nachlesen.

## Github Account einrichten
Um einen Github Account einzurichten, folgt man den Anweisungen [hier](https://github.com/join). Es gibt auch Alternativen dazu, aber Github wird am häufigsten
verwendet.

## Fork und Clone
Gehe auf die Website des Repositories des Numerikpraktikums und klicke auf den
Button `Fork` rechts oben. Dies erzeugt ein lokales Repository auf dem eigenen
Github Account. 

Im nächsten Schritt klont man das Repository auf den eigenen lokalen Rechner. 
Dazu klicke auf den grünen `Code` Button und kopiert die URL. Danach gibt man
den Befehl
```
git clone <URL des Repositories>
```
ein und wechselt mit
```
cd <repo-name>
```
in das lokale Verzeichnis.


## Arbeiten auf dem eigenen Branch
Es wird empfohlen immer einen eignen Branch für Code-Änderungen zu haben. So mit
kann der main-Branch möglichst fehlerfrei gehalten werden. Um dies zu machen,
nutz man den Befehl:
```
git checkout -b my-feature-branch
```
Dabei sollte man einen aussagekräftigen Namen wählen.

Nun kann man die gewünschten Änderungen einfügen in den Code. 

## Stage, Commit und Push
Bevor man einen Commit macht, müssen die Änderungen aufgenommen werden. Dazu
kann man einfach den Befehl
```
git add .
```
nutzen, der alle Dateien hinzufügt. Es ist aber auch möglich nur bestimmte
Dateien hinzuzufügen mit
```
git add <file-name>
```
Dabei können auch Wildcard-Expressions genutzt werden. 

Nun kann man einen einen Commit machen. Dies ist quasi ein Checkpoint für
Änderungen und an diesen Punkten, kann man die Versionen auch später noch
einsehen.
Dies geschieht über den Befehl
```
git commit -m "Commit Message"
```
Dabei ist eine Commit Message zwingend erforderlich.
Nun kann dies in das eigene Repository auf Github hochgeladen werden mit
```
git push origin my-feature-branch
```
Nun ist es im eigenen Repository auf Github. Wenn man nicht mit anderen Menschen
zusammenarbeitet, ist man nun fertig. Doch wir wollen ja zusammen arbeiten,
daher werden Pull-Requests gestellt. 

## Pull Request
Geh auf dein Repository auf Github und erstelle ein Pull Request. Dies erlaubt
den Eigentümer des ursprünglichen Repositories die Änderungen zu übernehmen.
Dabei muss man nur den Anweisungen auf der Website folgen.

## Up-to-Date bleiben
Diese Schritten können beliebig oft wiederholt werden. Allerdings kann sich
das Original-Repository auch verändern, während man selber lokal arbeitet. 
Daher muss es auch geupdated werden. Um eure lokaler Repository mit dem
originalen zu verknüpfen, benutzt den folgenden Befehl:
```
git remote add upstream <original-repo-url>
```

Dann kann man das lokale Repository updaten mit:
```
git fetch upstream
git checkout main
git merge upstream/main
```

Dies fügt die Verschiedenen Code-Versionen zusammen und es kann geupdatet werden
mit
```
git push origin main
```


