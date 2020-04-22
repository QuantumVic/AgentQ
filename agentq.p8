pico-8 cartridge // http://www.pico-8.com
version 20
__lua__
-- math
-- this code is part of qiskit.
--
-- copyright IBM 2020

-- custom math table for compatibility with the Pico8

math = {}
math.pi = 3.14159
math.max = max
math.sqrt = sqrt
math.floor = flr
function math.random()
  return rnd(1)
end
function math.cos(theta)
  return cos(theta/(2*math.pi))
end
function math.sin(theta)
  return -sin(theta/(2*math.pi))
end
function math.randomseed(time)
end
os = {}
function os.time()
end
-->8
-- microqiskit
-- this code is part of qiskit.
--
-- copyright ibm 2020

--math.randomseed(os.time())

function quantumcircuit ()

  local qc = {}

  local function set_registers (n,m)
    qc._n = n
    qc._m = m or 0
  end
  qc.set_registers = set_registers

  qc.data = {}

  function qc.initialize (ket)
    ket_copy = {}
    for j, amp in pairs(ket) do
      if type(amp)=="number" then
        ket_copy[j] = {amp, 0}
      else
        ket_copy[j] = {amp[0], amp[1]}
      end
    end
    qc.data = {{'init',ket_copy}}
  end

  function qc.add_circuit (qc2)
    qc._n = math.max(qc._n,qc2._n)
    qc._m = math.max(qc._m,qc2._m)
    for g, gate in pairs(qc2.data) do
      qc.data[#qc.data+1] = ( gate )    
    end
  end
      
  function qc.x (q)
    qc.data[#qc.data+1] = ( {'x',q} )
  end

  function qc.rx (theta,q)
    qc.data[#qc.data+1] = ( {'rx',theta,q} )
  end

  function qc.h (q)
    qc.data[#qc.data+1] = ( {'h',q} )
  end

  function qc.cx (s,t)
    qc.data[#qc.data+1] = ( {'cx',s,t} )
  end

  function qc.measure (q,b)
    qc.data[#qc.data+1] = ( {'m',q,b} )
  end

  function qc.rz (theta,q)
    qc.h(q)
    qc.rx(theta,q)
    qc.h(q)
  end

  function qc.ry (theta,q)
    qc.rx(math.pi/2,q)
    qc.rz(theta,q)
    qc.rx(-math.pi/2,q)
  end

  function qc.z (q)
    qc.rz(math.pi,q)
  end

  function qc.y (q)
    qc.z(q)
    qc.x(q)
  end

  return qc

end




function simulate (qc, get, shots)

  if not shots then
    shots = 1024
  end

  function as_bits (num,bits)
    -- returns num converted to a bitstring of length bits
    -- adapted from https://stackoverflow.com/a/9080080/1225661
    local bitstring = {}
    for index = bits, 1, -1 do
        b = num - math.floor(num/2)*2
        num = math.floor((num - b) / 2)
        bitstring[index] = b
    end
    return bitstring
  end

  function get_out (j)
    raw_out = as_bits(j-1,qc._n)
    out = ""
    for b=0,qc._m-1 do
      if output_map[b] then
        out = raw_out[qc._n-output_map[b]]..out
      end
    end
    return out
  end


  ket = {}
  for j=1,2^qc._n do
    ket[j] = {0,0}
  end
  ket[1] = {1,0}

  output_map = {}

  for g, gate in pairs(qc.data) do

    if gate[1]=='init' then

      for j, amp in pairs(gate[2]) do
          ket[j] = {amp[1], amp[2]}
      end

    elseif gate[1]=='m' then

      output_map[gate[3]] = gate[2]

    elseif gate[1]=="x" or gate[1]=="rx" or gate[1]=="h" then

      j = gate[#gate]

      for i0=0,2^j-1 do
        for i1=0,2^(qc._n-j-1)-1 do
          b1=i0+2^(j+1)*i1 + 1
          b2=b1+2^j

          e = {{ket[b1][1],ket[b1][2]},{ket[b2][1],ket[b2][2]}}

          if gate[1]=="x" then
            ket[b1] = e[2]
            ket[b2] = e[1]
          elseif gate[1]=="rx" then
            theta = gate[2]
            ket[b1][1] = e[1][1]*math.cos(theta/2)+e[2][2]*math.sin(theta/2)
            ket[b1][2] = e[1][2]*math.cos(theta/2)-e[2][1]*math.sin(theta/2)
            ket[b2][1] = e[2][1]*math.cos(theta/2)+e[1][2]*math.sin(theta/2)
            ket[b2][2] = e[2][2]*math.cos(theta/2)-e[1][1]*math.sin(theta/2)
          elseif gate[1]=="h" then
            for k=1,2 do
              ket[b1][k] = (e[1][k] + e[2][k])/math.sqrt(2)
              ket[b2][k] = (e[1][k] - e[2][k])/math.sqrt(2)
            end
          end

        end
      end

    elseif gate[1]=="cx" then

      s = gate[2]
      t = gate[3]

      if s>t then
        h = s
        l = t
      else
        h = t
        l = s
      end

      for i0=0,2^l-1 do
        for i1=0,2^(h-l-1)-1 do
          for i2=0,2^(qc._n-h-1)-1 do
            b1 = i0 + 2^(l+1)*i1 + 2^(h+1)*i2 + 2^s + 1
            b2 = b1 + 2^t
            e = {{ket[b1][1],ket[b1][2]},{ket[b2][1],ket[b2][2]}}
            ket[b1] = e[2]
            ket[b2] = e[1]
          end
        end
      end

    end

  end

  if get=="statevector" then
    return ket
  else

    probs = {}
    for j,amp in pairs(ket) do
      probs[j] = amp[1]^2 + amp[2]^2
    end

    if get=="fast counts" then

      c = {}
      for j,p in pairs(probs) do
        out = get_out(j)
        if c[out] then
          c[out] = c[out] + probs[j]*shots
        else
          if out then -- in case of pico8 weirdness
            c[out] = probs[j]*shots
          end
        end
      end
      return c

    else

      m = {}
      for s=1,shots do
        cumu = 0
        un = true
        r = math.random()
        for j,p in pairs(probs) do
          cumu = cumu + p
          if r<cumu and un then
            m[s] = get_out(j)
            un = false
          end
        end
      end

      if get=="memory" then
        return m

      elseif get=="counts" then
        c = {}
        for s=1,shots do
          if c[m[s]] then
            c[m[s]] = c[m[s]] + 1
          else
            if m[s] then -- in case of pico8 weirdness
              c[m[s]] = 1
            else
              if c["error"] then
                c["error"] = c["error"]+1
              else
                c["error"] = 1
              end
            end
          end
        end
        return c

      end

    end

  end

end
-->8
-- main
-- initialization

palt(0,false)
palt(4,true)

function _init()

	cartdata("agent_q_qiskit")
	
	--parameter init	
	iteration=0
	mapy=-48*8
	highscore=0

		
	init_parameters()

	music(1,100)
	lose_cursor_pos=1
	
	cls()
	title_screen()
	
	
	--quantum circuit init
	qc={}
	
	for num_qub=0,3 do
		setup_qc(num_qub,1)
	end
	
	--sprite init
	control_sprite=nil
	make_cursor()
	make_qubits()
	make_buttons()
	

	
	--game state init	
	intro_done=false
	game_started=false
	
	
	---------debug
	--game_started=true
	--intro_done = true
	--during_tutorial= false
	--gamemode=1
	--difficulty=1
	--lose_status=true
	---------end debug
end

function _update()
 scroll()
 
 --game state upd
	if script_active then
    script_update()
 end

	if not intro_done then
		script_run(intro)
		intro_done = true
	end
	
	if moving_tutorial then
		move_cursor()
	end
		
	--qiskit
	get_statevector()
	state_to_bool()
	
	if lose_status then
		lose_screen()
	end
	
		--game
	if game_started then
		move_cursor()
		game()
	 if timer_global>2 then
	 	shoot_qubit()
  end
  check_death()
  get_highscore()
	end
	

end

function _draw()
		cls()
  draw_map(mapy)
  
  --[[ debug draw
  print(qub[0].amp1)
  print(qub[0].amp2)
  print(qub[1].amp1)
  print(qub[1].amp2)
  print(qub[2].amp1)
  print(qub[2].amp2)
  --end debug draw]]

  if game_started then
  	draw_qubits(0)
			draw_qubits(1)
			draw_qubits(2)
			draw_cursor()
		end
  
		if game_started or lose_status then	
			draw_menu()
			draw_health()
			print("time",1,1)
			print(seconds,5,8)
			print("score",13,101)
			print(score,36,101)
			if difficulty==1 then
				print("high",13,110)
				print(highscore,36,110)
			end
		end
		
		  
  if lose_status then
  	draw_lose_screen()
  end
  
  if credits_status then
  	draw_credits()
  end
  
		if during_tutorial then
			draw_qubits(0)
			draw_qubits(1)
			draw_qubits(2)
			
			if control_sprite!=nil then
				control_show()
			end
			
			draw_cursor()
		end
		
    if text then
        rectfill(2,107,125,125,0)
        print(text, 3,108, text_color)
    end
    if responses then
        local top = 101 - 6 * #responses
        rectfill(70, top,
                 125, 105, 0)
        for i=1, #responses do
            print(responses[i],
                  72, top + i*6-4,
                  i==ans and 7 or 5)
        end
    end
    --screenshot()
end
-->8
-- map+qubit code


function draw_map(mapy)
	map(0,0,0,mapy,16,63)
end

function draw_lose_screen()
	for i=4,11 do
			spr(49,i*8,4*8)
			spr(52,i*8,9*8)
		for j=5,8 do
			spr(4,i*8,j*8)
		end
	end
	
	for j=5,8 do
			spr(54,3*8,j*8)
			spr(55,12*8,j*8)
	end
	spr(48,3*8,4*8)
	spr(50,12*8,4*8)
	spr(51,3*8,9*8)
	spr(53,12*8,9*8)
	
	print("agent q",64-7*2,4.5*8)
	
	print("play again",64-10*2,6.5*8+2)
	print("title screen",64-12*2,7.5*8+2)
	print("credits",64-7*2,8.5*8+2)
	
	draw_lose_cursor()
	
end

function make_qubits()
	qub={}
	qub[0]={x=16/2-1.5,y=2,size=3,sup=false,state=false,amp1=0,amp2=0}
	qub[1]={x=3,y=7,size=3,sup=false,state=false,amp1=0,amp2=0}
	qub[2]={x=16-6,y=7,size=3,sup=false,state=false,amp1=0,amp2=0}
	qub[3]={x=14,y=14,size=1,sup=false,state=false,amp1=0,amp2=0}
	
	state={}
	state[true]={sprite=17}
	state[false]={sprite=18}
end

function make_buttons()
	xbutton={sprite=90,size=3,x=15/2+4.5,y=12,changed=84}
	zbutton={sprite=93,size=3,x=15/2+2,y=12,changed=87}
end

function draw_qubits(q_num)
	for i=0,qub[q_num].size-1 do
		for j=0,qub[q_num].size-1 do
			spr(state[qub[q_num].state].sprite+bool_to_number[qub[q_num].sup]*2,
				(qub[q_num].x+i)*8,
				(qub[q_num].y+j)*8,
				1,1)
		end
	end
	spr(state[qub[q_num].state].sprite+16+bool_to_number[qub[q_num].sup]*2,
				(qub[q_num].x+1)*8,
				(qub[q_num].y+1)*8,
				1,1)
				
				
end

function draw_menu()
	for i = 10,13 do
			spr(4,i*8,13*8)
			spr(49,i*8,12*8)
			spr(52,i*8,14*8)
	end
	
	spr(48,9.5*8,12*8)
	spr(50,14*8,12*8)
	spr(51,9.5*8,14*8)
	spr(53,14*8,14*8)
	spr(54,9.5*8,13*8)
	spr(55,14*8,13*8)
	
	for i = 2,5 do
			spr(4,i*8,13*8)
			spr(49,i*8,12*8)
			spr(52,i*8,14*8)
	end
	
	spr(50,(15-9.5)*8,12*8)
	spr(48,(15-14)*8,12*8)
	spr(53,(15-9.5)*8,14*8)
	spr(51,(15-14)*8,14*8)
	spr(55,(15-9.5)*8,13*8)
	spr(54,(15-14)*8,13*8)
	
	spr(xbutton.sprite,xbutton.x*8,
	xbutton.y*8,xbutton.size,
	xbutton.size)
	spr(zbutton.sprite,zbutton.x*8,
	zbutton.y*8,zbutton.size,
	zbutton.size)
	
	if timer_global>2 then
		if btnp(5) and game_started then
			frame_temp_xbutton=timer_global
			xbutton.sprite=xbutton.changed
		end
		
		if frame_temp_xbutton!=nil 
			and timer_global==frame_temp_xbutton+3 then
				xbutton.sprite=90
		end	
		
		if btnp(4) then
			frame_temp_zbutton=timer_global
			zbutton.sprite=zbutton.changed	
		end
		
		if frame_temp_zbutton!=nil 
			and timer_global==frame_temp_zbutton+3 then
				zbutton.sprite=93
		end	
	end

end

function scroll()
	iteration+=1
	mapy=-48*8+iteration
	
	if (mapy==0) then
		mapy=-48*8
		iteration=0
	end
end

function draw_health()
	for i=1,health do
		spr(56,128-8*i,0)
	end
end
-->8
-- cursor code

function make_cursor()
	p={}
	p.qub=0
	p.sprite=5
	p.keys=0
	p.size=3
end

function draw_cursor()
		p.x=qub[p.qub].x
		p.y=qub[p.qub].y
		spr(p.sprite,p.x*8,
		p.y*8,p.size,p.size)
end

function draw_lose_cursor()
	spr(57,27,45+8*lose_cursor_pos)
end

function move_cursor()
		if btnp(⬆️) then
			p.qub=0
		elseif (btnp(⬅️) or (p.qub==0 and btnp(⬇️))) then
			p.qub=1
		elseif btnp(➡️) then
			p.qub=2
		end
end
-->8
-- dialog code

-- scripting variables
-------------------------------
text = nil
text_color = 7
responses = nil
ans = 1
routine = nil
script_active = false

-- initiate a script
function script_run(func)
    routine = cocreate(function()
        script_active = true
        func()
        script_active = false
    end)
    coresume(routine)
end

-- this is called every frame
-- and player input is ignored,
-- as long as there is a script
-- active.
function script_update()
    coresume(routine)
end


-- script commands
-------------------------------

function reveal_text(str)
    text = ""
    for i=1, #str do
        text = text..sub(str,i,i)
        yield()	
    end
end

function say(str)
    reveal_text(str)
    repeat
     yield()
    until btnp(5)
    text = nil
end

function say_z(str)
    reveal_text(str)
    repeat
        yield()
    until btnp(4)
    text = nil
end

function announce(str)
    text = str
    text_color = 12
    repeat
        yield()
    until btnp(5)
    text = nil
    text_color = 7
end

function ask(str, ...)
    reveal_text(str)
    responses = {...}
    ans = 1
    repeat
        yield()
        if btnp(2) and ans > 1 then
            ans -= 1
        elseif btnp(3) and ans < #responses then
            ans += 1
        end
    until btnp(5)
    text = nil
    responses = nil
end

-- execute multiple script
--  functions at once.
-- the main script resumes once
--  all functions are complete
function simultaneously(...)
    local routines = {}
    for f in all{...} do
        add(routines, cocreate(f))
    end
    repeat
        yield()
        local complete = true
        for c in all(routines) do
            if coresume(c) then
             complete = false
            end
        end
    until complete
end


function control_show()
		spr(8,qub[control_sprite].x*8,
						qub[control_sprite].y*8,
						3,3)
end

-- introduction script
intro = function ()
say [[hello, agent q.
how are you feeling today?]]
say [[good, i hope, because
today's mission is
extremely important.]]
say [[today, you are getting rid of
the errors inside our
quantum computer!]]
ask ([[do you know the details
of the mission, or do you wish
to go through the tutorial?]],"i'm good","please help!")
if ans==2 then
say[[ok, agent q.
pay close attention!]]
--qub[1].state = not qub[1].state
--qub[2].state = not qub[2].state
qc[0].x(0)
music(2)
during_tutorial=true
moving_tutorial=true
say[[on your screen, you have
three qubits of our
quantum computer.]]
say[[as you can see, the top qubit
is in the state |1>...]]
say[[...and there are two
more "empty" qubits, in
the state |0>.]]
say[[you can navigate through them
using the ⬆️,⬇️,⬅️,➡️ keys.
try it!]]
say[[let's use a cnot gate,
in order to "repeat" the
top qubit's state.]]
p.sprite=8
repeat
say[[i have given you a control,
please apply it to the
top qubit with ❎.]]
until (btnp(❎) and p.qub==0)
control_sprite=0
p.sprite=11
say[[now, here you have the other
part of the cnot gate:
the x gate.]]
repeat
say[[place the x gate on one
of the qubits to entangle
with ❎.]]
until (btnp(❎) and (p.qub==1 or p.qub==2))
control_sprite=nil
first_ent_qubit=p.qub
qc[p.qub].x(0)
p.sprite=5
say[[congratulations! you have
successfully used an 
entangling gate...]]
say[[...to repeat the state of 
the qubit.]]
say[[now, let's repeat the state
in the remaining qubit.]]
p.sprite=8
repeat
say[[please, place this control
on one of the qubits in
the state |1>.]]
until (btnp(❎) and (p.qub==0 or p.qub==first_ent_qubit))
control_sprite=p.qub
p.sprite=11
repeat					
say[[and now, please place
the x gate on the
remaining qubit.]]
until (btnp(❎) and (p.qub!=0 and p.qub!=first_ent_qubit))
qc[p.qub].x(0)
control_sprite=nil
p.sprite=5
say[[the interesting thing about
this procedure is that
it can be done...]]
say[[with any arbitrary
state of a qubit!]]
say[[now... why are we
doing this?]]
say[[from time to time, one of
the qubits may randomly
flip	its state!]]
say[[
this is called a "bit flip"]]
qc[0].x(0)
sfx(20)

say[[right there,
see what i told you?!]]
say[[that's a huge issue, and this
is where you come in!]]
say[[it is your duty to correct
these errors, agent q!]]
wait(2)
repeat
say[[select the flipped qubit
and press ❎ to flip it back
to its original state!]]
until (btnp(❎) and p.qub==0)
qc[p.qub].x(0)
sfx(18)

say[[
great!]]
say[[since we have three qubits,
if one of them suffers a bit
flip, we can compare it...]]
say[[...with the other two,
and correct the error!]]
say[[the best thing about this
kind of error correction is
that it works for any...]]
say[[...arbitrary quantum state.]]

say[[now, the qubits can be in
what we call 'superposition'.]]
say[[this means that they can be
in a state |+>=|0>+|1>...]]
say[[or also |->=|0>-|1>,
among an infinity of possible
states.]]

say[[let me quickly switch our
state to |+>.]]

qc[0].x(0)
qc[1].x(0)
qc[2].x(0)
hadamard()

say[[the same way that a qubit
can randomly flip its
state...]]
say[[...the sign of the phase may
also randomly change.]]
say[[
this is called a
phase flip.]]

z_gate(1)
sfx(20)

say[[damn it!
there it is!]]

repeat

say_z[[quick! select the flipped
sign qubit and press z to
correct its phase!]]
until (btnp(4) and p.qub==1)

z_gate(p.qub)
sfx(18)

say[[
phew! that was close...]]
say[[the good news is,
i think you're getting
the gist of it!]]
say[[it's time for you to
go out there and correct all
those pesky quantum errors.]]
say[[
but be careful!]]
say[[one error may me manageable,
but if you stack up two errors
at the same time...]]
say[[it will be impossible for us
to tell right from wrong and
come back to the right state.]]
say[[oh! and every time you apply
the wrong gate, or to the wrong
qubit, you will lose hearts.]]
say[[also, keep in mind the quantum
computer will be used while you
correct the errors...]]
say[[so you may see all the qubits
change together. just keep
correcting the odd one out!]]
say[[it is your duty to keep
our quantum state alive.
we trust in you, agent q.]]
say[[
don't let us down!]]
moving_tutorial=false
end
during_tutorial=false
	qc={}
	
	for num_qub=0,3 do
		setup_qc(num_qub,1)
		qc[num_qub].x(0)
	end
ask([[choose the difficulty.]],
"normal","easy")

difficulty = ans

say[[
good luck then, agent q.]]

wait(1)
game_started=true
music(-1)
music(0)
wait(14)
end

-->8
-- game code

function title_screen()
	highscore=dget(0)
	repeat
	
		print("▥▥▥▥▥▥▥▥",64-16*2,61-18)
		print("▥            ▥",64-16*2,61-24)
		print("▥  agent  q  ▥",64-16*2,61-30)
		print("▥            ▥",64-16*2,61-36)
		print("▥▥▥▥▥▥▥▥",64-16*2,61-42)
		
		print("♥ vic pina ♥",64-14*2,61-6)
		print("press ❎",64-8*2,61+40)
		print("high score   ",64-13*2,61+32)

		print(highscore,64+10*2,61+32)
		
	until btnp(❎)
	
end

function lose_screen()

	lose_cursor()
	timer_global=0
	
	if btnp(❎) and lose_cursor_pos==1 then
		--play again
		restart_game()
	elseif btnp(❎) and lose_cursor_pos==2 then
	 --title screen
		run()
	elseif btnp(❎) and lose_cursor_pos==3 then
		credits_status=true 
	end
	
end

function draw_credits()

	for i=2,13 do
			spr(49,i*8,3*8)
			spr(52,i*8,14*8)
		for j=4,13 do
			spr(4,i*8,j*8)
		end
	end
	
	for j=4,13 do
			spr(54,1*8,j*8)
			spr(55,14*8,j*8)
	end
	spr(48,1*8,3*8)
	spr(50,14*8,3*8)
	spr(51,1*8,14*8)
	spr(53,14*8,14*8)
	
	print("game....vicente pina",64-20*2,32)
	print("music...vicente pina",64-20*2,42)
	
	print("..acknowledgements..",64-20*2,64)
	print("james wootton",64-13*2,82)
	print("huang junye",64-11*2,90)
	print("samanvay sharma",64-15*2,98)
	print("geckojsc",64-8*2,106)
	
	print("press 🅾️",64-8*2,14)
	
	if btnp(4) then
		credits_status=false
	end
end

function restart_game()
		wait(1)
		init_parameters()
		
		game_started=true
		lose_status=false
		
		for num_qub=0,3 do
		setup_qc(num_qub,1)
		end	
	
		x_all()
		music(0)
end

function lose_cursor()
	if btnp(⬆️) then
		lose_cursor_pos-=1
	elseif btnp(⬇️) then
		lose_cursor_pos+=1
	end
	
	if lose_cursor_pos==0 then
		lose_cursor_pos=3
	elseif lose_cursor_pos==4 then
		lose_cursor_pos=1
	end	
end

function bit_flip(prob,bit_min,bit_period)
if timer_bit > bit_min  and timer_bit%bit_period==0 and math.random()<prob then
		if math.random()<0.33 then
			qc[0].x(0)
			sfx(20)
		elseif math.random()<0.67 then
			qc[1].x(0)
			sfx(20)
		else
			qc[2].x(0)
			sfx(20)
		end
	end
	timer_bit+=1
end

function phase_flip(prob,phase_min,phase_period)
if timer_phase > phase_min  and timer_phase%phase_period==0 and math.random()<prob then
		if math.random()<0.33 then
			z_gate(0)
			sfx(20)
		elseif math.random()<0.67 then
			z_gate(1)
			sfx(20)
		else
			z_gate(2)
			sfx(20)
		end
	end
	timer_phase+=1
end

function init_parameters()
	timer_bit=0
	timer_phase=0
	timer_global=0
	seconds=0
	score=0
	health=3
	flip_period=150
end



function measure_qc()
	 qc.measure(0,0)
 	qc.measure(1,1)
 	qc.measure(2,2)
end

function shoot_qubit()
--bit flip if not in
--superposition and x press
--phase flip if in
--superposition and z press

if btnp(5) then
	if qub[3].sup then
		health-=1
		sfx(17)
		xbutton.changed=132
	elseif qub[3].state==qub[p.qub].state then
		health-=1
		sfx(17)
		xbutton.changed=132
	else
		qc[p.qub].x(0)
		xbutton.changed=84
		sfx(18)
		score+=1
	end
end

if btnp(4) then
	if not qub[3].sup then
		health-=1
		sfx(17)
		zbutton.changed=135
	elseif qub[3].state==qub[p.qub].state then
		health-=1
		sfx(17)
		zbutton.changed=135
	else
		z_gate(p.qub)
		zbutton.changed=87
		sfx(18)
		score+=1
	end
end
end


function game()
	if not qub[3].sup then
		bit_flip(1,1,flip_period)	
	else
		phase_flip(1,1,flip_period)
	end
	
	if difficulty==1 then
		decrease_period()
	end
	
 hadamard_flip(1,150,200)
 x_flip(1,100,205)
 z_flip(1,175,210)
	timer_global+=1
	seconds=math.floor(timer_global/30)
end

function decrease_period()
	if seconds>3 and flip_period>40 then
		if timer_global%30==0 then
		 flip_period-=1
		end
	end
end

function hadamard_flip(prob,h_min,h_period) 
	if timer_global%h_period==0
	and math.random()<prob
	and timer_global>h_min then
		hadamard()
		sfx(19)
	end
end


function x_flip(prob,x_min,x_period) 
	if timer_global%x_period==0
	and math.random()<prob
	and timer_global>x_min 
	and not qub[3].sup then
		x_all()
		sfx(19)
	end
end


function z_flip(prob,z_min,z_period) 
	if timer_global%z_period==0
	and math.random()<prob
	and timer_global>z_min 
	and qub[3].sup then
		z_all()
		sfx(19)
	end
end

function check_death()
	error_count=0
	
	if qub[3].state != qub[0].state then
		error_count+=1
	end
	if qub[3].state != qub[1].state then
		error_count+=1
	end
	if qub[3].state != qub[2].state then
		error_count+=1
	end
	
	if error_count>=2 
	or health==0 then
	lose_game()
	end
end

function get_highscore()
	if difficulty==1
	and score > highscore then
		highscore = score
		dset(0,highscore)
	end
end

function lose_game()
	game_started=false
	lose_status=true
	xbutton.changed=90
	zbutton.changed=93
	music(2)
end

bool_to_number={ [true]=1, [false]=0 }

function wait(a) for i = 1,a do flip() end end

function screenshot()
	cls()
			draw_map(mapy)
  	draw_qubits(0)
			draw_qubits(1)
			draw_qubits(2)
			draw_cursor()
			
			
			for i=4,11 do
			spr(49,i*8,12*8-4)
			spr(52,i*8,14*8-4)
			spr(4,i*8,13*8-4)
	end
	
			spr(54,3*8,13*8-4)
			spr(55,12*8,13*8-4)
	spr(48,3*8,12*8-4)
	spr(50,12*8,12*8-4)
	spr(51,3*8,14*8-4)
	spr(53,12*8,14*8-4)
	
	print("▥  agent q  ▥",64-15*2,13*8-3)

end
-->8
--qiskit code


function setup_qc(qc_num,num_qub)
	qc[qc_num] = quantumcircuit()
	qc[qc_num].set_registers(num_qub)
end

function z_gate(num_qub)
	qc[num_qub].h(0)
	qc[num_qub].x(0)
	qc[num_qub].h(0)
end

function hadamard()
	qc[0].h(0)
	qc[1].h(0)
	qc[2].h(0)
	qc[3].h(0)
end


function x_all()
	qc[0].x(0)
	qc[1].x(0)
	qc[2].x(0)
	qc[3].x(0)
end


function z_all()
	qc[0].h(0)
	qc[0].x(0)
	qc[0].h(0)
	qc[1].h(0)
	qc[1].x(0)
	qc[1].h(0)
	qc[2].h(0)
	qc[2].x(0)
	qc[2].h(0)
	qc[3].h(0)
	qc[3].x(0)
	qc[3].h(0)
end

--statevector simulation
function get_statevector()
	result ={}
	for num_qub=0,3 do
		result[num_qub] = simulate(qc[num_qub],"statevector")
		qub[num_qub].amp1=result[num_qub][1][1]
		qub[num_qub].amp2=result[num_qub][2][1]
	end
end

function state_to_bool()
	for i=0,3 do
		if qub[i].amp1 == 1 then
			qub[i].state=false
			qub[i].sup=false
		elseif qub[i].amp1==0 then
			qub[i].state=true
			qub[i].sup=false
		elseif qub[i].amp2>0.6 and qub[i].amp2<0.8 then
			qub[i].state=false
			qub[i].sup=true
		elseif qub[i].amp2<-0.6 and qub[i].amp2>-0.8 then
			qub[i].state=true
			qub[i].sup=true
		end			
	end
end
__gfx__
00000000777777773333333344444444000000001111111144444444111111114444444444444444444444444444444444444444444444440000000000000000
00000000777777773333333344444444000000001777766144444444166777714444444444444444444444444444444444444444444444440000000000000000
00700700777777773333333344444444000000001777766144444444166777714444444444444444444444444400000000000000000000440000000000000000
0007700077777777333333334444444400000000177111114444444411111771444444444444444444444444440cccccccccccccccccc0440000000000000000
0007700077777777333333334444444400000000177144444444444444441771444444444444444444444444440cccccccccccccccccc0440000000000000000
0070070077777777333333334444444400000000166144444444444444441661444444444000000444444444440cccccccccccccccccc0440000000000000000
00000000777777773333333344444444000000001661444444444444444416614444444000cccc0004444444440ccc11cccccccc11ccc0440000000000000000
000000007777777733333333444444440000000011114444444444444444111144444400cccccccc00444444440ccc111cccccc111ccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee4444444444444444444444444444440cccccccccc0444444440cccc111cccc111cccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee4444444444444444444444444444400cccccccccc0044444440ccccc111cc111ccccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee444444444444444444444444444440cccccccccccc044444440cccccc111111cccccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee444444444444444444444444444440cccccccccccc044444440ccccccc1111ccccccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee444444444444444444444444444440cccccccccccc044444440ccccccc1111ccccccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee444444444444444444444444444440cccccccccccc044444440cccccc111111cccccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee4444444444444444444444444444400cccccccccc0044444440ccccc111cc111ccccc0440000000000000000
00000000ccccccccaaaaaaaabbbbbbbbeeeeeeee4444444444444444444444444444440cccccccccc0444444440cccc111cccc111cccc0440000000000000000
00000000ccc11cccaa9999aabbbbbbbbeeeeeeee11114444444444444444111144444400cccccccc00444444440ccc111cccccc111ccc0440000000000000000
00000000cc111ccca999999abbbbbbbbeee88eee1661444444444444444416614444444000cccc0004444444440ccc11cccccccc11ccc0440000000000000000
00000000cc111ccca99aa99abbbbbbbbeee88eee166144444444444444441661444444444000000444444444440cccccccccccccccccc0440000000000000000
00000000ccc11ccca99aa99ab333333be888888e177144444444444444441771444444444444444444444444440cccccccccccccccccc0440000000000000000
00000000ccc11ccca99aa99ab333333be888888e177111114444444411111771444444444444444444444444440cccccccccccccccccc0440000000000000000
00000000ccc11ccca99aa99abbbbbbbbeee88eee1777766144444444166777714444444444444444444444444400000000000000000000440000000000000000
00000000c111111ca999999abbbbbbbbeee88eee1777766144444444166777714444444444444444444444444444444444444444444444440000000000000000
00000000c111111caa9999aabbbbbbbbeeeeeeee1111111144444444111111114444444444444444444444444444444444444444444444440000000000000000
33333333333333333333333330000000000000000000000330000000000000034444444444444444000000000000000000000000000000000000000000000000
30000000000000000000000330000000000000000000000330000000000000034004004444477444000000000000000000000000000000000000000000000000
30000000000000000000000330000000000000000000000330000000000000030880880444477744000000000000000000000000000000000000000000000000
30000000000000000000000330000000000000000000000330000000000000030888880444477774000000000000000000000000000000000000000000000000
30000000000000000000000330000000000000000000000330000000000000030888880444477774000000000000000000000000000000000000000000000000
30000000000000000000000330000000000000000000000330000000000000034088804444477744000000000000000000000000000000000000000000000000
30000000000000000000000330000000000000000000000330000000000000034408044444477444000000000000000000000000000000000000000000000000
30000000000000000000000333333333333333333333333330000000000000034440444444444444000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00333300000030000033330000033300000330000033300000033000003303000000300000300300003030000033330000333300000303000003300000030300
00300000000030000003000000300300003003000003030000300300000030000033330000333300003333000030030000300300000303000030030000330300
00030000000003000033330000300300003003000003030000030000003030000030030000300300003003000030000000333300000333000030000000033300
00030000003003000003000000300300003330000033300000300000003030000003000000030000003003000003000000300300000300000003000000030300
00003000030000300003000000300300003000000003000000300300000030000003000000030000000003000003030000300300000300000000300000330300
00003000033333300000330000033300000333000003000000033000003303000030330000003300000330000000330000303000003333000033330000030300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00033000003000000000030000000000444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00300300003000000000300000300300444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00300300003000000000030000300300444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00300300003000000000300000030300444433333333333333334444444433333333333333334444444466666666666666664444444466666666666666664444
0030030000030000000300000000330044443bbbbbbbbbbbbbb3444444443bbbbbbbbbbbbbb34444444467777777777777764444444467777777777777764444
0003300000003300003333000033030044443bbbbbbbbbbbbbb3444444443bbbbbbbbbbbbbb34444444467777777777777764444444467777777777777764444
0000000000000000000000000000000044443bb33bbbbbb33bb3444444443bbb33333333bbb34444444467766777777667764444444467776666666677764444
0000000000000000000000000000000044443bb333bbbb333bb3444444443bb3333333333bb34444444467766677776667764444444467766666666667764444
0000000000000000000000000000000044443bbb333bb333bbb3444444443bbbbbbbb333bbb34444444467776667766677764444444467777777766677764444
0000000000000000000000000000000044443bbbb333333bbbb3444444443bbbbbbb333bbbb34444444467777666666777764444444467777777666777764444
0000000000000000000000000000000044443bbbbb3333bbbbb3444444443bbbbbb333bbbbb34444444467777766667777764444444467777776667777764444
0000000000000000000000000000000044443bbbbb3333bbbbb3444444443bbbbb333bbbbbb34444444467777766667777764444444467777766677777764444
0000000000000000000000000000000044443bbbb333333bbbb3444444443bbbb333bbbbbbb34444444467777666666777764444444467777666777777764444
0000000000000000000000000000000044443bbb333bb333bbb3444444443bbb333bbbbbbbb34444444467776667766677764444444467776667777777764444
0000000000000000000000000000000044443bb333bbbb333bb3444444443bb333333333bbb34444444467766677776667764444444467766666666677764444
0000000000000000000000000000000044443bb33bbbbbb33bb3444444443bbb333333333bb34444444467766777777667764444444467776666666667764444
0000000000000000000000000000000044443bbbbbbbbbbbbbb3444444443bbbbbbbbbbbbbb34444444467777777777777764444444467777777777777764444
0000000000000000000000000000000044443bbbbbbbbbbbbbb3444444443bbbbbbbbbbbbbb34444444467777777777777764444444467777777777777764444
00000000000000000000000000000000444433333333333333334444444433333333333333334444444466666666666666664444444466666666666666664444
00000000000000000000000000000000444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00000000000000000000000000000000444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00000000000000000000000000000000444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
00000000000000000000000000000000444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
c400b400346400001400005400008400444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
00000000000000000000000000000000444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
b400d400006400000400006400000000444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
00000000000000000000000000000000444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
0000e400000015003400007400000000444488888888888888884444444488888888888888884444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448eeeeeeeeeeeeee8444444448eeeeeeeeeeeeee84444000000000000000000000000000000000000000000000000
d40035f400000500240000840054000044448eeeeeeeeeeeeee8444444448eeeeeeeeeeeeee84444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448ee88eeeeee88ee8444444448eee88888888eee84444000000000000000000000000000000000000000000000000
b40015d40000000000d40000007400d444448ee888eeee888ee8444444448ee8888888888ee84444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448eee888ee888eee8444444448eeeeeeee888eee84444000000000000000000000000000000000000000000000000
940000b40000350000c40000009400e444448eeee888888eeee8444444448eeeeeee888eeee84444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448eeeee8888eeeee8444444448eeeeee888eeeee84444000000000000000000000000000000000000000000000000
840000f40000250000b4e40000b400f444448eeeee8888eeeee8444444448eeeee888eeeeee84444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448eeee888888eeee8444444448eeee888eeeeeee84444000000000000000000000000000000000000000000000000
34000005a4001500005414000000000544448eee888ee888eee8444444448eee888eeeeeeee84444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448ee888eeee888ee8444444448ee888888888eee84444000000000000000000000000000000000000000000000000
2400000094000500740424000000001544448ee88eeeeee88ee8444444448eee888888888ee84444000000000000000000000000000000000000000000000000
0000000000000000000000000000000044448eeeeeeeeeeeeee8444444448eeeeeeeeeeeeee84444000000000000000000000000000000000000000000000000
140000008400f400640034005400002544448eeeeeeeeeeeeee8444444448eeeeeeeeeeeeee84444000000000000000000000000000000000000000000000000
00000000000000000000000000000000444488888888888888884444444488888888888888884444000000000000000000000000000000000000000000000000
6400a40094000000c400440064000035444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
00000000000000000000000000000000444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
00008400a40000002500000074000000444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
00000000000000000000000000000000444444444444444444444444444444444444444444444444000000000000000000000000000000000000000000000000
00005400b40000150000000084000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000034000000000500840000d4004400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000024e400e400f400440035f4006400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000000d400a400000034002515007494000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000400c40084000000140005350000a4000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00333300000030000033330000033300000330000033300000033000003303000000300000300300003030000033330000333300000303000003300000030300
001400b400740000000400e4000000b4000000000000000000000000000000000033330000333300003333000030030000300300000303000030030000330300
00030000000003000033330000300300003003000003030000030000003030000030030000300300003003000030000000333300000333000030000000033300
00240000009400b4000000c400840000000000000000000000000000000000000000000000000000003003000003000000300300000300000003000000030300
00003000030000300003000000300300003000000003000000300300000030000003000000030000000003000003030000300300000300000000300000330300
0034005400a400c40000000000b40034000000000000000000000000000000000030330000003300000330000000330000303000003333000033330000030300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0044006400b400d40000050000e40054000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00033000003000000000030000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000007400c400e40000f40025050074000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00300300003000000000030000000300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000084000000050000c40035150064000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00300300000300000003000000003300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00d400a4000000f4d400940014350044000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00f4009400000000a400640074250024000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00e400000000e4007400340064000004000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000d4000000040054000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000064000025f40000c4000044000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000840000350500001400003400b400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000a40000041500000400000000f400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000c4000014000000c4000000007400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000e400004400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__label__
00300300000000000030000000030300000000000000000000000000000000000000000000030300000000000000000000000000000030000000000000030300
00300000000000000030000000033300000000000000000000000000000000000000000000033300000000000000000000000000003030000000000000033300
00030000000000000030000000030000000000000000000000000000000000000000000000030000000000000000000000000000003030000000000000030000
00030300000000000003000000030000000000000000000000000000000000000000000000030000000000000000000000000000000030000000000000030000
00003300000000000000330000333300000000000000000000000000000000000000000000333300000000000000000000000000003303000000000000333300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00300300000000000000000000333300000000000000000000000000000000000000000000333300000000000000000000000000003003000000000000033000
00333300000000000000000000300300000000000000000000300300000000000000000000300300000000000000000000000000003333000000000000300300
00300300000000000000000000300000000000000000000000300300000000000000000000333300000000000000000000000000003003000000000000300000
00030000000000000000000000030000000000000000000000030300000000000000000000300300000000000000000000000000000300000000000000030000
00030000000000000000000000030300000000000000000000003300000000000000000000300300000000000000000000000000000300000000000000003000
00003300000000000000000000003300000000000000000000330300000000000000000000303000000000000000000000000000000033000000000000333300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00003000000000000000000000030300000000000000000000000300000000000000000000333300000330000000000000000000003333000000000000030300
003333000000000000000000003303000000000000000000000011111111aaaaaaaa111111110300003003000000000000000000003003000000000000330300
003003000000000000000000000333000000000000000000000017777661aaaaaaaa166777710000003000000000000000000000003000000000000000033300
000300000000000000000000000303000000000000000000000017777661aaaaaaaa166777710000000300000000000000000000000300000000000000030300
000300000000000000000000003303000000000000000000000317711111aaaaaaaa111117710300000030000000000000000000000303000000000000330300
00303300000000000000000000030300000000000000000000331771aaaaaaaaaaaaaaaa17713300003333000000000000000000000033000000000000030300
00000000000000000000000000000000000000000000000000001661aaaaaaaaaaaaaaaa16610000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000001661aaaaaaaaaaaaaaaa16610000000000000000000000000000000000000000000000000000
00033300000000000000000000033000003030000000000000301111aaaaaaaaaaaaaaaa11113000000030000000000000000000000000000000000000033000
0030030000000000000000000030030000333300000000000030aaaaaaaaaa9999aaaaaaaaaa0300000030000000000000000000000000000000000000300300
0030030000000000000000000030030000300300000000000030aaaaaaaaa999999aaaaaaaaa0300000003000000000000000000000000000000000000300300
0030030000000000000000000030030000300300000000000030aaaaaaaaa99aa99aaaaaaaaa3000003003000000000000000000000000000000000000300300
0030030000000000000000000030030000000300000000000003aaaaaaaaa99aa99aaaaaaaaa0000030000300000000000000000000000000000000000300300
0003330000000000000000000003300000033000000000000000aaaaaaaaa99aa99aaaaaaaaa0000033333300000000000000000000000000000000000033000
0000000000000000000000000000000000000000000000000000aaaaaaaaa99aa99aaaaaaaaa0000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000aaaaaaaaa999999aaaaaaaaa0000000000000000000000000000000000000000000000000000
0033330000000000000000000000000000300300000000000003aaaaaaaaaa9999aaaaaaaaaa3300003333000000000000000000000000000000000000300000
00030000000000000000000000000000003333000000000000301111aaaaaaaaaaaaaaaa11110000000300000000000000000000000000000000000000300000
00333300000000000000000000000000003003000000000000301661aaaaaaaaaaaaaaaa16610000003333000000000000000000000000000000000000300000
00030000000000000000000000000000000300000000000000301661aaaaaaaaaaaaaaaa16610000000300000000000000000000000000000000000000300000
00030000000000000000000000000000000300000000000000301771aaaaaaaaaaaaaaaa17713000000300000000000000000000000000000000000000030000
000033000000000000000000000000000000330000000000000317711111aaaaaaaa111117713000000033000000000000000000000000000000000000003300
000000000000000000000000000000000000000000000000000017777661aaaaaaaa166777710000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000017777661aaaaaaaa166777710000000000000000000000000000000000000000000000000000
000030000000000000000000000000000000300000000000000311111111aaaaaaaa111111110000000333000000000000333000000000000000000000000300
00003000000000000000000000000000003333000000000000330300000000000030030000000000003003000000000000030300000000000000000000003000
00000300000000000000000000000000003003000000000000033300000000000003000000000000003003000000000000030300000000000000000000000300
00300300000000000000000000000000000300000000000000030300000000000030000000000000003003000000000000333000000000000000000000003000
03000030000000000000000000000000000300000000000000330300000000000030030000000000003003000000000000030000000000000000000000030000
03333330000000000000000000000000003033000000000000030300000000000003300000000000000333000000000000030000000000000000000000333300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00033000000000000030300000000000003003000000000000000000000000000033330000000000000330000000000000033000000000000000000000000000
00300300000000000033330000000000003333000000000000000000000000000030030000000000003003000000000000300300000000000000000000300300
00030000000000000030030000000000003003000000000000000000000000000033330000000000003003000000000000030000000000000000000000300300
00300000000000000030030000000000000300000000000000000000000000000030030000000000003330000000000000300000000000000000000000030300
00300300000000000000030000000000000300000000000000000000000000000030030000000000003000000000000000300300000000000000000000003300
00033000000000000003300000000000000033000000000000000000000000000030300000000000000333000000000000033000000000000000000000330300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000300000000000003030000000000000000000000000000000030000000000000000000000000000330300000000000000000000000000
000000000000000000333300cccccccccccccccccccccccc00000000000000000000300000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000300300cccccccccccccccccccccccc00000000000000000000030000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000030000cccccccccccccccccccccccc00000000000000000000300000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000030000cccccccccccccccccccccccc00000000000000000003000000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000303300cccccccccccccccccccccccc00000000000000000033330000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000000000cccccccccccccccccccccccc00000000000000000000000000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000000000cccccccccccccccccccccccc00000000000000000000000000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000333000cccccccccccccccccccccccc00000000003000000000000000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000030300ccccccccccc11ccccccccccc00000000003000000000000000000000ccccccccccc11ccccccccccc000000000000000000000000
000000000000000000030300cccccccccc111ccccccccccc00000000003000000000000000000000cccccccccc111ccccccccccc000000000000000000000000
000000000000000000333000cccccccccc111ccccccccccc00000000003000000000000000000000cccccccccc111ccccccccccc000000000000000000000000
000000000000000000030000ccccccccccc11ccccccccccc00000000000300000000000000000000ccccccccccc11ccccccccccc000000000000000000000000
000000000000000000030000ccccccccccc11ccccccccccc00000000000033000000000000000000ccccccccccc11ccccccccccc000000000000000000000000
000000000000000000000000ccccccccccc11ccccccccccc00000000000000000000000000000000ccccccccccc11ccccccccccc000000000000000000000000
000000000000000000000000ccccccccc111111ccccccccc00000000000000000000000000000000ccccccccc111111ccccccccc000000000000000000000000
000000000000000000033300ccccccccc111111ccccccccc00000000000330000000000000003000ccccccccc111111ccccccccc000000000003300000000000
000000000000000000300300cccccccccccccccccccccccc00000000003003000000000000333300cccccccccccccccccccccccc000000000030030000000000
000000000000000000300300cccccccccccccccccccccccc00000000003003000000000000300300cccccccccccccccccccccccc000000000030030000000000
000000000000000000300300cccccccccccccccccccccccc00000000003003000000000000030000cccccccccccccccccccccccc000000000033300000000000
000000000000000000300300cccccccccccccccccccccccc00000000003003000000000000030000cccccccccccccccccccccccc000000000030000000000000
000000000000000000033300cccccccccccccccccccccccc00000000000330000000000000303300cccccccccccccccccccccccc000000000003330000000000
000000000000000000000000cccccccccccccccccccccccc00000000000000000000000000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000000000cccccccccccccccccccccccc00000000000000000000000000000000cccccccccccccccccccccccc000000000000000000000000
000000000000000000333300cccccccccccccccccccccccc00000000000303000000000000033000cccccccccccccccccccccccc000000000003300000000000
00000000000000000003000000300300000000000030030000000000003303000000000000300300000000000030030000330300000000000030030000000000
00000000000000000033330000300000000000000030000000000000000333000000000000300300000000000030030000033300000000000003000000000000
00000000000000000003000000030000000000000003000000000000000303000000000000333000000000000003030000030300000000000030000000000000
00000000000000000003000000003000000000000000300000000000003303000000000000300000000000000000330000330300000000000030030000000000
00000000000000000000330000333300000000000033330000000000000303000000000000033300000000000033030000030300000000000003300000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000030300000000000030300000000000000000000000000000033300000000000000030000300000000000000033030000300300
00000000000000000000000000030300000000000033330000000000000000000000000000300300000000000000300000300000000000000000300000333300
00000000000000000000000000033300000000000030030000000000000000000000000000300300000000000000030000300000000000000030300000300300
00000000000000000000000000030000000000000030030000000000000000000000000000300300000000000000300000300000000000000030300000030000
00000000000000000000000000030000000000000000030000000000000000000000000000300300000000000003000000030000000000000000300000030000
00000000000000000000000033333333333333333333333333333333333333333333333333333333333333333333333333333333000000000033030000003300
00000000000000000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000000000
00000000000000000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000000000
00000000003333000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000303000
00000000003000000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000333300
00000000000300000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000300300
00000000000300000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000300300
00000000000030000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000000300
00000000000030000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000033000
00000000000000000000000030000000007070707000000000777007707770770077700000070000000000707070700000000003000000000000000000000000
00000000000000000000000030000000007070707000000000707070007000707007000000707000000000707070700000000003000000000000000000000000
00000000000030000000000030000000007070707000000000777070007700707007000000707000000000707070700000000003000000000000000000333300
00000000000030000000000030000000007070707000000000707070707000707007000000770000000000707070700000000003000000000000000000300300
00000000000003000000000030000000007070707000000000707077707770707007000000077000000000707070700000000003000000000000000000300000
00000000003003000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000030000
00000000030000300000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000030300
00000000033333300000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000003300
00000000000000000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000000000
00000000000000000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000000000000000000000000
00000000003333000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000030000000000000000000
00000000000300000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003003333000000000000000000
00000000003333000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003003003000000000000000000
00000000000300000000000030000000000000000000000000000000000000000000000000000000000000000000000000000003000300000000000000000000
00000000000300000000000033333333333333333333333333333333333333333333333333333333333333333333333333333333000300000000000000000000
00000000000033000000000000000000000000000000330000000000000033000000000000000000000000000030300000000000003033000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000333000000000000333000000000000030300000000000003333000000000000000000000000000000000000000000003333000000000000033300
00000000003003000000000000030300000000000033330000000000003003000000000000000000000000000000000000000000003003000000000000300300
00000000003003000000000000030300000000000030030000000000003333000000000000000000000000000000000000000000003000000000000000300300
00000000003003000000000000333000000000000030030000000000003003000000000000000000000000000000000000000000000300000000000000300300
00000000003003000000000000030000000000000000030000000000003003000000000000000000000000000000000000000000000303000000000000300300
00000000000333000000000000030000000000000003300000000000003030000000000000000000000000000000000000000000000033000000000000033300
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000330000000000000033000000000000033330000000000000303000000000000000000000330000000000000000000000330000000000000333000

__gff__
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
d81401800000000000008080487880487420690a000000000000800a6e6f75740a0a6e6f000000000000807270650a692232350a788048808048786d7370096100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__map__
0040004c00480000004100505300004a00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0041004b004700000040004e0000004b00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004200000049004b0000004c0048000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00430045004a004c00000000004b004300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00440046004b004d00005000004e004500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000047004c004e00004f005250004700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000000480000005000004c005351004600000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004d004a0000004f4d0049004153004400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004f0049000000004a0046004752004200000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004e000000004e00470043004600004000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000004d00000040004500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000460000524f00004c00004400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00004800005350000041000043004b0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00004a00004051000040000000004f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00004c0000410000004c00000000470000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0000004c004500500046004b4300004600000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0050004b4700004f0045004a4200004900000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004f000048000000004400004000004a00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004e000049004600000000004100004b00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
004d0000000044004600000000004b4c00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
474c00005000430047004f0000494a0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
484b00004f00000049004e000048410000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
49000000524e000048004d4d0047000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000004d00004a00004c0042440000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00530000004c00490000004b0043420000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00520047004b00480043004a0000430000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00510048004a0047004400000000004f00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
0050444a00000000004500004700004e00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
53004500000052000046000048004b4d00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
52004900004751000000004449004a0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
4e004a0044450000420000430000490000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__sfx__
011700001c0331c0331f6351c0331c03317053170531c63510053100531365510053100531705317053106551b0331b0331e6351b0331b03317053170531b6350f0530f053126550f0530f05317053170530f655
01170000107351373517735137351073513735177351373510735137351773513735107351373517735137350f7351273517735127350f7351273517735127350f7351273517735127350f735127351773512735
01170000102351023510003102351360510235102351023510003100031000310205136051020510205102350f2350f235100030f235136050f2350f2350f2351000310003100030f205136050c0030c0030f335
011d00001c0231c0231f6251c0231c02317023170231c62510023100231362510023100231702317023106251b0231b0231e6251b0231b02317023170231b6250f0230f023126250f0230f02317023170230f625
011d0000107351373517735137351073513735177351373510735137351773513735107351373517735137350f7351273517735127350f7351273517735127350f7351273517735127350f735127351773512735
011d00001021510215100031021513605102151021510215000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
011700000000000000000000000000000000000000000000233241f3241f324233241f324233241f324000000000000000000000000000000000000000000000232241e2241e224232241e224232241e22400000
011d00000000000000000000000000000000000000000000233141f3141f314233141f314233141f314000000000000000000000000000000000000000000000232141e2141e214232141e214232141e21400000
011d00000000000000000000000000000000000000000000233141f3141f314233141f314233141f3140000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
011000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
01170000182530c7000c5000c4000c3000c2002410024005240032460324503244032430324103240032420324203242031820300000000000000000000000000000000000000000000000000000000000000000
011000002405524005000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
010500002415126153000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
010400002435126553000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
010600001030010200101001000010300104001050000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
__music__
03 00010206
03 03044507
03 03040508

