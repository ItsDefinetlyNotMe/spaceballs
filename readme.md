After installing the necessary libraries with:
```
pip install -r requirements.txt
```
The program can be started in two modes, either in the editor with:
```
Python main.py editor /scene/example.scene
```
which opens an editor window which lets you adjust the scene while plotting:
![image of the editor](./images/editor.png)

Or in plot mode, which simply plots the scene:
```
Python main.py plot /scene/example.scene 
```
![image of a plot](./images/example_figure.png)
in either options more details can be printed with the 'print details' button or the '-p' option respectively.

Furthermore the option -e allows for the setting of an epsilon that must be between collisions, this can be set to 0 to allow the ball to collide twice or more with the same planet at once.
This can also be set to regulate float imprecision, which might allow for otherwise illegal collisions 