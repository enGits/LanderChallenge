

## File Structure:
```
/game/
│
├── index.html
├── main.js
│
├── /scenes/
│   ├── menu.js
│   ├── levelMenu.js
│   ├── tutorial.js
│   ├── gameOver.js
│   └── orbitGame.js
│
├── /assets/
│   └── logo.png
|
├── /core/
│   ├── body.js
│   ├── spacecraft.js
│   ├── surfaceFeature.js
│   ├── mountain.js
```


```
// base
class SurfaceFeature
class Mountain extends SurfaceFeature

// physical bodies
class Body
class Spacecraft extends Body
```
