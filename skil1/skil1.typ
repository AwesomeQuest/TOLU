#figure(image("ru-logo.svg", width: 40%))

#align(center, text(20pt)[#smallcaps[
  Háskólinn í Reykjavík
]])

#align(center, text(17pt)[#smallcaps[
  Skilaverkenfi
]])

#align(center, pad(y:6em,text(26pt)[
  *Skilaverkefni 1*
]),)

#pad(y:6em,grid(
  columns: (1fr),
  align(center)[
    Torfi Þogrímsson \
    #link("mailto:torfi22@ru.is")
  ],
  
))


#grid(
  columns: (1fr),
  align(center)[
    Supervised by  \
    Ingunn Gunnarsdóttir \
    #link("mailto:ingunngunnars@ru.is") \ \
    Olivier Moschetta \
    #link("mailto:olivierm@ru.is")
  ],
  
)

#pagebreak()

#set text(lang:"is")

#set math.equation(numbering: none)

#let numeq(eq) = math.equation(block: true, numbering: "(1)", eq)

= Lýsing 

#numeq[$
  cos(k L) cosh(k L) = -1
$]<eq:main>

$
  y(x,t) = A(x) B(t) \
  A(x) = cosh(k x) - cos(k x) + (cos(k L) + cosh(k L))/(sin(k L) + sinh(k L)) (sin(k x) - sinh(k x)) \
  B(t) = delta/A(L) cos(omega t)
$


= Verkefni


Við gefum okkur að lengd  stangarinnar er $L = 30 "mm"$, massaþéttleiki stangarinnar er $λ = 0.8 "g"/"mm"^3$ og stífnistuðull er $E I = 1.09 times 1010 "Pa" dot "mm"^2$


 1. Plottið graf fallsins $f(x) = cos(x) cosh(x) + 1$ á bilinu [0, 10] svo að þrjár rætur sjáist vel. Við viljum finna rót f skv. jöfnunni (1)

#figure(image("Graphics/2024-11-25-12-20-35.png", width: 50%), caption: [Hér má sjá @eq:main á $sinh^(-1)$ scala.])

 2.  Notið helmingunaraðferð til að finna minnstu jákvæðu rót fallsins $f$ með 4 réttum aukastöfum. Reiknið í kjölfarið tíðni $ω_1$ sem samsvarar þessu gildi

```julia 
function bisect(f, a, b, tol)
    m = (a + b)/2
    N = ceil(Int, log2((b-a)/2tol))
    for i in 1:N
        a, b = f(a)*f(m) < 0 ? (a, m) : (m, b)
        m = (a + b)/2
    end
    m
end

bisect(x->cos(x)*cosh(x)+1, 0,2, 1e-5)
1.8750991821289062
```

$
  k L = 1.8750 & => k = 0.062503 "mm"^(-1)\
  k^4 = (omega^2 lambda)/(E I) &=> omega_1 = k^2 sqrt((E I)/(lambda)) = 0.1262 ["rad"/"s"]
$

#linebreak()

3.  Hve margar ítranir af helmingunarðaferð nægja til að fá 4 rétta aukastafi í lið 2? Sýnið fræðilega útreikninga og sannreynið með teljara inní lykkjunn

$
  N = 17
$

4. @eq:main hefur í raun óendanlega margar lausnir. Notið aðferð Newtons til að finna næst minnstu jákvæðu rót $f$ með 4 réttum aukastöfum. Rökstyðjið val á upphafsgildinu. Reiknið í kjölfarið næst minnstu sveiflutíðni $ω_2$

```julia
function newton(f, x0, tol)
    df(x) = ForwardDiff.derivative(f,x)
    prev = x0
    curr = x0 - f(x0)/df(x0)
    
    while abs(prev-curr) > tol
        prev,curr = curr, curr - f(curr)/df(curr)
    end
    curr
end

let
	EI = 1.09*1010Pa*mm^2
	λ  = 0.8g/mm^3
	k  = (newton(x->cos(x)*cosh(x)+1, 4.8, 1e-5) ± 5e-5)/30mm
	k^2*sqrt(EI/λ) |> Hz
end
0.908217 ± 1.9e-5 Hz

```
$
  omega_2 = 0.9082 ["rad"/"s"]
$

5.  Notið aðferð að eigin vali til að reikna tuttugu minnstu sveiflutíðnir stangarinnar. Hér þarf bæði að ganga úr skugga um að engin lausn hafi gleymst eða verið tvítali

Fyrst $cosh(x)$ er stranglega vaxandi fyrir öll $x>0$, þá er eina leiðin fyrir jöfnu (@eq:main) að hafa fleyri rætur $cos(x)$ partin að þvínga það í núll. $cos$ hefur lotuna $2pi$ og krossar $x$ ásinn með bilinu $pi$. Þar af leiðandi skoðum clamped (@eq:main) á bilinu $[0,20pi]$

#grid(align:center + horizon, columns: 2,
  [#figure(
    image("Graphics/2024-11-25-13-10-01.png"),
    caption: [Clamp-uð @eq:main]
  )<fig:clamp>],
  [#figure(
    image("Graphics/2024-11-25-13-26-50.png"),
    caption: [Afleiða Clamp-uð @eq:main],
    )<fig:difclamp>]

)


Í @fig:clamp og @fig:difclamp sjáum við allar 20 fyrstu ræturnar. 

Við notum þetta til að ákvarða $x_0$ fyrir `newton`

Sem milli skref skoðum við afleiðu fallsins @fig:difclamp, þar sem er búið að finna hápunkta og lágpunkta.
#figure(grid(align: center + horizon, columns: 2, 
pad(y:1em)[Rætur $k L$],
pad(y:1em)[Tíðnir $omega$ [Hz]],
pad(x:1em)[#table(align: right,  columns: 5,
..("1.8751","4.6941","7.8548","10.9955","14.1372","17.2788","20.4204","23.5619","26.7035","29.8451","32.9867","36.1283","39.2699","42.4115","45.5531","48.6947","51.8363","54.9779","58.1195","61.2611",).map(i=>[#i]))],
pad(x:1em)[#table(align: right,  columns: 5,
..("0.145","0.908","2.543","4.983","8.238","12.306","17.188","22.883","29.392","36.714","44.850","53.800","63.563","74.140","85.531","97.735","110.753","124.584","139.229","154.688").map(i=>[#i]))],
), caption: [Töflur af rótum á @eq:main $k L$ og samsvarandi tíðnum $omega$])

6.  Við notum gildi á $ω$ og $k$ sem hafa verið reiknuð í lið 4. Búið til hreyfimynd af stönginni. Hér má endilega hafa skjölin *hreyfi.m* og *hreyfi.py* á Canvas sem fyrirmynd

#figure(grid(align: center + horizon, columns: 3,
  ..("Graphics/2024-11-25-17-36-43.png",
      "Graphics/2024-11-25-17-38-58.png",
      "Graphics/2024-11-25-17-39-09.png",
      "Graphics/2024-11-25-17-39-18.png",
      "Graphics/2024-11-25-17-39-26.png",
      "Graphics/2024-11-25-17-39-35.png")
      .enumerate().map(((ind,value))=>pad(y:1em)[
        #figure(numbering: none,image(value), 
        caption: $t=display(T*ind/6)$)]) 
), caption: [Hér sést hliðrun stangarinnar í mismunandi tíma punktum, þar sem $T = display((2pi)/omega_2) approx 6.9182" s"$ er lota $B(t)$]
)