# III-EIIMFE

Repository for the EIIMFE

To run properly this repository we need the library [SpecialFucntions](https://juliamath.github.io/SpecialFunctions.jl/latest/). To install this you can type in your Julia terminal the symbol `]` to open the Pkg enviroment and type

`pkg> add SpecialFunctions`

in case of faulire visit the original webpage.

## The notebooks
I recomend to use the _Pluto notebooks_ included in this repositoy. [Pluto](https://www.juliapackages.com/p/pluto) is an interesting alternative to the [_Jupyter notebooks_](https://jupyter.org/) given that this is build specially for the Julia Language and is easy to manipulate and inderstand[[1](https://www.youtube.com/watch?v=IAF8DjrQSSk&t=5s)]. To easy install this package just type

`pkg> add Pluto`

into the pkg enviroment. And additionally install the `PlutoUI` typing into the same enviroment

`pkg> add PlutoUI`

Some of the examples need the `Plots` package, then for an interactive experiance install this package too

`pkg> add Plots`

Finally, to run the notebooks type into the julia termina

`julia> import Pluto`
`julia> Pluto.run()`

and your web browser will display some menu like this

![image](https://user-images.githubusercontent.com/60673156/142267004-12a413df-bd03-44cb-a4ab-e20a6a5d3e87.png)

you can open the notebooks just selecting this from the path were is allocated.
