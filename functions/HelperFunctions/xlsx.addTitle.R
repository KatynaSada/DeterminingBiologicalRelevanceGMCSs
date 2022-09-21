#http://www.sthda.com/english/wiki/r-xlsx-package-a-quick-start-guide-to-manipulate-excel-files-in-r

#++++++++++++++++++++++++
# Helper function to add titles
#++++++++++++++++++++++++
# - sheet : sheet object to contain the title
# - rowIndex : numeric value indicating the row to 
#contain the title
# - title : the text to use as title
# - titleStyle : style object to use for title
xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-createCell(rows, colIndex=1)
  setCellValue(sheetTitle[[1,1]], title)
  setCellStyle(sheetTitle[[1,1]], titleStyle)
}