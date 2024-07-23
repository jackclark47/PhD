####TODO:
###ADD: pause if spacebar is hit DONE ###
###ADD: control description at game start DONE ####
###ADD: game starts when Enter is pressed DONE ####
###ADD: game saves and quits when backspace is pressed DONE ####
###ADD: game loads the last save when an option in play() is selected. DONE ####
###ADD: a greeting that asks a player if they'd like to start a new game or continue a save. DONE ####

library(dplyr)
library(grid)
library(tcltk)


#Apple constructor
new_Apple <- function(height, width){
  apple_loc <- c(sample(1:height, 1), sample(1:width, 1))
  structure(list(coords = apple_loc) ,class = 'Apple')
}

#Defining the Game class
new_Game <- function(height, width){
  if (!is.numeric(height) | !is.numeric(width)) stop('height and width must be numeric')
  board <- structure(list(height = height, width = width), class = 'Board')
  
  init_body <- list(
    c(8,1),
    c(7,1), 
    c(7,2), 
    c(7,3),
    c(6,3),
    c(5,3)
  )
  
  snake <- new_Snake(init_body = init_body, init_direction = c(0,1))
  while(TRUE){
    apple <- new_Apple(height, width)
    if(!all(apple$coords %in% init_body)){
      break
    }
  }
  
  game <- structure(list(
    board = board, 
    snake = snake,
    apple = apple,
    
    DISPLAY_CHARS = list(
      EMPTY = "#003366",
      BODY = "#00cc00",
      HEAD = '#ffff66',
      APPLE = '#ff0000'
    ),
    
    DIR_RIGHT = c(0,1),
    DIR_LEFT = c(0,-1),
    DIR_UP = c(-1,0),
    DIR_DOWN = c(1,0),
    
    INPUT_CHAR_UP = '87',
    INPUT_CHAR_DOWN = '83',
    INPUT_CHAR_LEFT = '65',
    INPUT_CHAR_RIGHT = '68',
    
    SCORE = 0
  ),
  class = 'Game')
  return(game)
}

#Produce a matrix representing the board
game_matrix <- function(game){
  board <- game$board
  snake <- game$snake
  matrix <- matrix(data = game$DISPLAY_CHARS$EMPTY, nrow = board$height, ncol = board$width)
  for(part in snake$body){
    matrix[part[1], part[2]] <- game$DISPLAY_CHARS$BODY
  }
  head = head(snake)
  matrix[head[[1]][1], head[[1]][2]] <- game$DISPLAY_CHARS$HEAD
  matrix[game$apple[[1]][1], game$apple[[1]][2]] <- game$DISPLAY_CHARS$APPLE
  return(matrix)
}

#Render the matrix to the console
game_render <- function(game){
  matrix <- game_matrix(game)
  grid.raster(matrix, interpolate = FALSE)
}

#Obtain the next position for the head to move in to
next_position <- function(position, step){
  return(
    c(position[[1]] + step)
  )
}

#Constructor for the Snake class
new_Snake <- function(init_body, init_direction){
  structure(list(body = init_body, direction = init_direction), class = 'Snake')
}

#Makes the next step of the snake
take_step <- function(game, position){
  game$snake$body <- game$snake$body[-1]
  game$snake$body <- append(game$snake$body, list(c(position[1], position[2])))
  return(game)
}

#Sets the direction for the next movement of the snake
set_direction <- function(game, direction){
  game$snake$direction <- direction
  return(game)
}

#gets location of the snakes head
head <- function(snake){
  return(snake$body[length(snake$body)])
}

#grows the snake once an apple is eaten
grow <- function(game, position){
  game$snake$body <- append(game$snake$body, list(c(position[1], position[2])))
  return(game)
}

#Function to save the current game state
save_game <- function(game){
  h <- paste('Height: ', game$board$height, '\n', sep = '')
  w <- paste('Width: ', game$board$width, '\n', sep = '')
  sc <- paste('Score: ', game$SCORE, '\n', sep = '')
  ap <- paste('Apple: ', game$apple$coords[1], ' ', game$apple$coords[2], '\n', sep = '')
  d <- paste('Direction: ', game$snake$direction[1], ' ', game$snake$direction[2], '\n', sep = '')
  b <- paste('Body:', game$snake$body, '\n', sep = '', collapse = '')
  writable<- paste(h, w, sc, ap, d, b, sep = '')
  writable <- substr(writable, 1, nchar(writable)-1)
  write(writable, 'D:/snakesave.txt')
}

load_game <- function(){
  data <- readLines('D:/snakesave.txt')
  height <- as.numeric(gsub("\\D", "", data[1]))
  width <- as.numeric(gsub("\\D", "", data[2]))
  score <- as.numeric(gsub("\\D", "", data[3]))
  apple <- gsub("\\D", "", data[4])
  apple <- c(substr(apple, 1, 1), substr(apple, 2,2)) %>%
    as.numeric()
  
  direction <- gsub("\\D", "", data[5])
  direction <- c(substr(direction,1,1), substr(direction, 2,2)) %>%
    as.numeric()
  
  body <- list(gsub("\\D", "", data[6:length(data)]))
  bodylist <- list()
  for(i in 1:length(body[[1]])){
    print(i)
    vec <- c(substr(body[[1]][i], 1,1), substr(body[[1]][i], 2,2)) %>%
      as.numeric() %>%
      list()
    bodylist[i] <- vec
    }
  body = bodylist

  game <- new_Game(height = height, width = width)
  game$SCORE <- score
  game$snake$direction <- direction
  game$snake$body <- body
  game$apple$coords <- apple
  return(game)
}

play <- function(game, edge_collision = 'OFF'){
  
  message('Welcome to SNAKE')
  continue <- readline(prompt = 'Press 1 to start a new game. Press 2 to continue from a saved game')
  if(continue == '2'){
    print(game)
    game <- load_game()
    print(game)
  }
  message('Controls:')
  message('Movement: wasd')
  message('Quit: Escape')
  message('Pause: Space')
  message('Save and quit: Backspace')
  invisible(readline(prompt = "Press Return to start snake.\n"))
  message("Game started.")
  
  if (.Device == "null device") dev.new()
  tt <- tcltk::tktoplevel()
  invisible(tcltk::tkwm.title(tt, "Snake control"))
  invisible(tcltk::tkwm.geometry(tt, "300x100+9000+500")) #https://stackoverflow.com/questions/14910858/how-to-specify-where-a-tkinter-window-opens
  SetFocus(tt)
  game_render(game)
  input <- '65'
  

  pause = FALSE
  print('========')
  print(game)
  while(TRUE){
    invisible(tkbind(tt, '<Key>', function(k) {input <<- k}))
    print(input)
    if(is.na(input)){
      Sys.sleep(1)
      next
    }
    
    if(input == '8'){
      print(game)
      save_game(game)
      message('Game saved. Quitting...')
      break
    }
    
    
    if(input == '32'){
      input = NA
      if (pause) message("Pause released.") else message("Game paused.")
      flush.console()
      print(pause)
      pause <- !pause
    }
    if(!pause){
      if(!is.na(input)){
        if(input == game$INPUT_CHAR_UP && !all(game$snake$direction == game$DIR_DOWN)){
          game <- set_direction(game, game$DIR_UP)
          
        }
        if(input == game$INPUT_CHAR_DOWN && !all(game$snake$direction == game$DIR_UP)){
          game <- set_direction(game, game$DIR_DOWN)
          
        }
        if(input == game$INPUT_CHAR_LEFT && !all(game$snake$direction == game$DIR_RIGHT)){
          game <- set_direction(game, game$DIR_LEFT)
        }
        if(input == game$INPUT_CHAR_RIGHT && !all(game$snake$direction == game$DIR_LEFT)){
          game <- set_direction(game, game$DIR_RIGHT)
        }
      }

      
      next_position <- next_position(head(game$snake), game$snake$direction)
      
      if(edge_collision == 'ON'){
        if(0 %in% next_position | (game$board$height+1) %in% next_position[1] | (game$board$width+1) %in% next_position[2]){
          break
        }
      } else{
        
        if(0 %in% next_position[1]){
          next_position[1] <- game$board$height
        }
        
        if(0 %in% next_position[2]){
          next_position[2] <- game$board$width
        }
        
        if((game$board$height+1) %in% next_position[1]){
          next_position[1] <- 1
        }
        
        if((game$board$width+1) %in% next_position[2]){
          next_position[2] <- 1
        }
      }
      
      if(list(next_position) %in% game$snake$body){
        break
      }
      
      if(all(next_position %in% game$apple$coords)){
        game$SCORE <- game$SCORE + 1
        
        game <- grow(game, next_position)
        
        while(TRUE){
          game$apple <- new_Apple(game$board$height, game$board$width)
          if(!all(game$apple$coords %in% game$snake$body)){
            break
          }
        }
      } else{
        game <- take_step(game, next_position)
      }
    }
    
    
    
    game_render(game)
    Sys.sleep(0.1)
  }
  print('Game over!')
  print(paste("You scored: ", game$SCORE, ' points.', sep = ''))
}

#Make new Game object and render it
game <- new_Game(20,20)
play(game, edge_collision = 'OFF')
