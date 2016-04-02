# Graphics_MIS_PT
CIS560 project - Multiple Importance Sampling (MIS) path tracer.


## Sample Renders

### Multiple Importance Sampling 
Different light source types:

|Sphere|![](Rendering/balanced_veach_direct.bmp)|
|-----|-------|
|Cube|![](Rendering/cube.bmp)|
|Ring|![](Rendering/ring.bmp)|
<br>

### Specular Transmission BxDF Implementation

|Perfect Reflection|
|------------------|
|![](Rendering/transmission.bmp)|
<br>

### Transmission and specular weighted material Implementation

|Glass1|Glass2|
|------|------|
|![](Rendering/glass.jpg)|![](Rendering/specular_transmission.jpg)|
<br>

### Perfect reflection bxdf implemantation<br>

|Perfect Reflection|
|------------------|
|![](Rendering/specular.bmp)|
<br>

### Direct lighting sample VS MIS<br>

|Direct Lighting Only|MIS|
|------|------|
|![](Rendering/directManySpheres.bmp)|![](Rendering/indirectManySpheres.bmp)|
