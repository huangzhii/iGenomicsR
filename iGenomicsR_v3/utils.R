smartModal <- function(error=c(T,F), title = "Title", content = "Content"){
  if(error){
    showModal(modalDialog(
      title = title, footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p(content),
          style = "margin: auto; text-align: center"
      )
    ))
  }
  else{
    showModal(modalDialog(
      title = title, footer = NULL,
      div(class = "busy",
          p(content),
          img(src="images/loading.gif"),
          style = "margin: auto; text-align: center"
      )
    ))
  }
}