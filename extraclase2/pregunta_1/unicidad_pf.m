

function c = unicidad_pf(f, a, b)
         syms x
         f_sym = sym(f);
         f_d = diff(f_sym,x);
         d = existencia_pf(f_sym, f_d , a, b); #Paso 1: demostrar existencia
         f_dd = diff(f_d, x); #Paso 2: 2da derivada de F
         if !d
           printf('La funcion no tiene existencia \n')
           c = 0;
         else
           printf('La funcion tiene existencia \n')
           f_d_num = function_handle(f_d);
           puntos_criticos = solve( f_dd == 0, x); #Paso 3: Encontar puntos criticos
           puntos_reales = [];

            for i = 1:length(puntos_criticos)
                punto = double(puntos_criticos(i));
                if isreal(punto) && imag(punto) == 0 && punto >= a && punto <= b # Filtro de reales y dentro del intervalo [a, b]
                    puntos_reales = [puntos_reales, real(punto)];
                endif
            endfor

         endif
         puntos_evaluados = [];
         #Paso 4: Evaluar los puntos en f'

         puntos_evaluados = [puntos_evaluados, f_d_num(puntos_reales), f_d_num(a),  f_d_num(b)];
         #Paso 5: Asignar Maximo y Minimo y verificar que pertenecen a ]-1,1[
         ydmax = max(puntos_evaluados);
         ydmin = min(puntos_evaluados);

         if -1 < ydmin && ydmin < 1 && -1 < ydmax && ydmax < 1
           printf('La funcion tiene un punto fijo \n')
           c = 1
         else
           c = 0
         endif












function d = existencia_pf(f_sym, f_d, a, b)
          #Paso 1 Derivar F, que ya lo hicimos
          syms x
          f = function_handle(f_sym);
          try

            puntos_criticos = solve( f_d == 0, x); #Paso 2 Encontar puntos criticos

            puntos_reales = [];

            for i = 1:length(puntos_criticos) # Paso 3 Encontar puntos criticos
                punto = double(puntos_criticos(i));
                if isreal(punto) && imag(punto) == 0 && punto >= a && punto <= b # Filtro de reales y dentro del intervalo [a, b]
                    puntos_reales = [puntos_reales, real(punto)];
                endif
            endfor

            puntos_evaluados = [];
            puntos_evaluados =  [puntos_evaluados, f(puntos_reales)];
            if  a<= f(a) && f(a) <= b
              puntos_evaluados = [puntos_evaluados, f(a)];
            endif

             if  a<= f(b) && f(b) <= b
              puntos_evaluados = [puntos_evaluados, f(b)];
            endif

            ymin = min( puntos_evaluados); #Paso 4 Asignar los ymax y ymin y verificar
            ymax = max( puntos_evaluados);


            if  a<= ymin && ymin <= b  && a<= ymax && ymax <= b #Verificacion final
              d = 1;
            else
              d = 0;
            endif

          catch err
            d = 0;
          end_try_catch









