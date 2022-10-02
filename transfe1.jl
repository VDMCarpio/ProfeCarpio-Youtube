### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 22a63557-2f3a-4818-a814-e3de952d33e0
md"""
# Ejemplo de transferencia de calor en dos dimensiones por conducción en estado transitorio usando diferencias finitas."""

# ╔═╡ 4794a12f-dc6b-455b-aa3c-20f000d28040
md"""
## Descripción del problema

Considere un sistema en dos dimensiones comprendido por un área cuadrada con una longitud de $3$ cm, donde las condiciones de frontera son:

**Superior (en $y=L$):** Temperatura constante de 100°C.

**Inferior (en $y=0$):** Temperatura constante de 300°C.

**Derecha (en $x=L$):** Convección con un coeficiente convectivo de $45 \frac{W}{m^2 K}$ y una temperatura de fluido de 20°.

**Izquierda (en $x=0$):** Calor constante de $1000 \frac{W}{m^2}$
"""

# ╔═╡ f5ff7b10-3e6a-11ed-11c8-7788c0ec33ce
md"""
## CONDICIONES DE FRONTERA

En este bloque se definen las condiciones a las cuales se encuentra expuesto el sistema analizado.
"""

# ╔═╡ 40736bb8-c58e-428d-a5c0-fbb216f003b1
begin
	# Condición inicial
	Tin = 20. #°C
	
	# Condición de frontera superior (y=alto)
	Ts = 100. #°C
	
	# Condición de frontera inferior (y=0)
	Ti = 300. #°C
	
	# Condición de frontera derecha (x=ancho)
	h = 45. #W/m^2-K
	Tamb = 20. #°C
	
	# Condición de frontera izquierda (x=0)
	q = 1000. #W/m^2

	# Longitud en x
	Lx = 3/100 #m

	# Longitud en y. Como se considera un sistema cuadrado se pude definir igual a la longitud en x
	Ly = Lx
end

# ╔═╡ c44adec1-eba8-46f3-b88d-0de0e89a51ee
md"""
## Propiedades

Aquí se definen las propiedades del material que serán usadas en el problema.
"""

# ╔═╡ 71d8ff3f-cb26-4993-8b19-e6aa7d0c1cf0
begin
	# Conductividad térmica
	k = 20 #W/mK

	# Difusividad térmica
	α = 6.694*10^(-6) #m^2/s
end

# ╔═╡ 6e9e767d-1e5c-4007-97a7-5b899e5e0cfa
md"## Condiciones numéricas

En esta sección se definen las condiciones de la solución numérica."

# ╔═╡ 2e39463c-cc41-4822-a421-cc1539cd9bf1
begin
	# Se define el número de nodos en la dirección de x y y.
	nx = 5
	ny = nx # Se decide usar la misma cantidad de nodos que en x.
	
	# Se evalúa la distancia entre nodos. Para este ejemplo consideramos que la distancia en x será igual a las distancias en y.
	Δx = Lx / (nx - 1)
	Δy = Δx 

	# Definimos el diferencial del tiempo. En este caso se comprobó en clase que este diferencial mantiene el criterio de estabilidad para la solución explícita.
	Δt = 2 #s

	# Evaluamos los números adimensionales ya discretizados.
	Fo = α*Δt/(Δx*Δy)
	Bi = h*Δx/k

	# Finalmente definimos el tiempo final del cálculo.
	tf = 10 # Tiempo en segundos.
	
	# El número de pasos en el tiempo, para efectos de usar la nomenclatura usada en clase, lo definimos como p.
	p = trunc(Int,tf/Δt) # Aquí estoy usando la función trunc, esto sólo es para asegurar que la variable p sea Entera, ya que como son pasos no tendremos valores racionales.
end

# ╔═╡ 722270e1-fe51-40bc-96c5-6de647d2b068
md"""
## Inicializando vectores y matrices
"""

# ╔═╡ 08288960-d4a6-487b-aebc-523fc0940d04
md"""
Inicializando la matriz donde se almacenarán los valores del tienpo viejo.
"""

# ╔═╡ 998ec869-bc56-4630-8347-e39f947b96da
Tvm=fill(Tin,(nx,ny))

# ╔═╡ 4eb666b6-ce3c-4009-89d0-f3cbd3ad9f3c
md"""
Inicializando la matriz donde se almacenarán los valores de la temperatura del tiempo nuevo.
"""

# ╔═╡ ac69c46f-c6e5-4f67-9361-98296bf313fd
Tnm=fill(0.,(nx,ny)) # Inicializando la matriz de temperatura en el tiempo nuevo.

# ╔═╡ e329675a-ed82-43b1-ab21-fe9407e49f4f
md"""
En la siguiente parte del código se evaluarán las ecuaciones en cada paso del tiempo, donde estas ecuaciones son:

**Nodos intermedios**

$T^{p+1}_{i,j} = Fo \left( T^{p}_{i+1,j} + T^{p}_{i-1,j} + T^{p}_{i,j-1} + T^{p}_{i,j+1}\right) + \left( 1-4Fo\right) T^{p}_{i,j}$

**Condición de frontera superior** la cual es de temperatura constante $Ts$

$T^{p+1}_{1,j} = Ts$

**Condición de frontera inferior** que es de temperatura constante $Ts$

$T^{p+1}_{ny,j} = Ti$

**Condición de frontera de la izquierda** la cual es la condición de calor constante

$q'' = -k \frac{ΔT}{Δx}$

que despejada para la temperatura de los nodos de la frontera izquierda queda como

$T^{p}_{i,1} = \frac{q'' Δx}{k} + T^{p}_{i,2}$

Si consideramos un nodo fantasma a la izquierda ($T_{i,0}$) de la superficie la ecuación se puede expresar como

$q'' = -k \frac{T^{p}_{i,2} - T^{p}_{i,0}}{2Δx}$

Si sustituimos en la ecuación de nodos intermedios (cuando suponemos que se evalúa en los nodos de la superficie izquierda, donde $T_{i,1}$) se obtiene la ecuación de la condición de frontera.

$T^{p+1}_{i,j} = Fo \left( T^{p}_{i+1,j} + T^{p}_{i-1,j} + 2T^{p}_{i,j-1}\right) + \left( 1-4Fo\right) T^{p}_{i,j} + 2\frac{Fo Δx q''}{k}$

**Condición de frontera de la derecha** la cual es por convección, de acuerdo a la tabla obtenida del libro del incrópera, para esta condición de frontera podemos usar la siguiente expresión:

$T^{p+1}_{i,nx} = Fo\left(2T^{p}_{i,nx-1} + T^{p}_{i+1,nx} + T^{p}_{i-1,nx} +2 Bi T_{∞}\right) + \left(1-4Fo - 2 Bi Fo\right) T^{p}_{i,nx}$

NOTA: El subíndice $i$ indica los renglones y el $j$ las columnas acorde al orden de `Julia`.

"""

# ╔═╡ 15285887-f17d-4247-be45-f1f3187d411e
begin
	for time in 1:p
		
		i = ny
		for j in 1:nx
			# Aplicando condición de frontera superior
			Tnm[1,j] = Ts
			# Aplicando condición de frontera inferior
			Tnm[ny,j] = Ti
		end
		local j = 1
		for i in 2:ny-1
			# Aplicando condición de la izquierda (flujo de calor constante)
			#Tnm[i,j] = q*Δx/k + Tnm[i,j+1]
			Tnm[i,j] = Fo*(Tvm[i+1,j]+Tvm[i-1,j]+2*Tvm[i,j+1]) + (1-4*Fo)*Tvm[i,j] + 2*Fo*Δx*q/k
		end
		j = nx
		for i in 2:ny-1
			# Aplicando condición de la derecha (convección)
			Tnm[i,j] = Fo*(2*Tvm[i,j-1] + Tvm[i+1,j] + Tvm[i-1,j] + 2*Bi*Tamb) + (1 - 4*Fo - 2*Bi*Fo)*Tvm[i,j]
		end
		for i in 2:ny-1
			for j in 2:nx-1
				Tnm[i,j]=Fo*(Tvm[i,j+1]+Tvm[i,j-1]+Tvm[i+1,j+1]+Tvm[i-1,j]) + (1-4*Fo)*Tvm[i,j]
			end
		end

		if time != p Tvm = Tnm end #Reiniciando la solución a tiempo viejo
	end
	Tnm
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─22a63557-2f3a-4818-a814-e3de952d33e0
# ╟─4794a12f-dc6b-455b-aa3c-20f000d28040
# ╟─f5ff7b10-3e6a-11ed-11c8-7788c0ec33ce
# ╠═40736bb8-c58e-428d-a5c0-fbb216f003b1
# ╟─c44adec1-eba8-46f3-b88d-0de0e89a51ee
# ╠═71d8ff3f-cb26-4993-8b19-e6aa7d0c1cf0
# ╟─6e9e767d-1e5c-4007-97a7-5b899e5e0cfa
# ╠═2e39463c-cc41-4822-a421-cc1539cd9bf1
# ╟─722270e1-fe51-40bc-96c5-6de647d2b068
# ╟─08288960-d4a6-487b-aebc-523fc0940d04
# ╠═998ec869-bc56-4630-8347-e39f947b96da
# ╟─4eb666b6-ce3c-4009-89d0-f3cbd3ad9f3c
# ╠═ac69c46f-c6e5-4f67-9361-98296bf313fd
# ╟─e329675a-ed82-43b1-ab21-fe9407e49f4f
# ╠═15285887-f17d-4247-be45-f1f3187d411e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
