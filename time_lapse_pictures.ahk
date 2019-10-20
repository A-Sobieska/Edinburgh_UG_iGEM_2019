; Taking pictures with a laptop web camera for a time lapse for the monitoring of hydrogen collection
; Zeke Rice-Gray
; for University of Edinburgh Undergraduate iGEM team 2019
; AutoHotKey script designed for Windows 10

#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
#Warn  ; Enable warnings to assist with detecting common errors.
#SingleInstance ; Do not run multiple instances at the same time
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.

SetTitleMatchMode 3 ; Exact matching mode for WinActivate

Loop
{
	Sleep (10 * 60 * 1000) ; sleep for 10 minutes
	WinActivate Camera ; choose the window with the camera
	Sleep 100 ; small delay
	Send {Space} ; take picture
}
