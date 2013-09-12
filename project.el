;; see ~/.emacs.d/emacs-rc-cedet.el
;; this file should be included to ~/.emacs.d/emacs-rc-cedet.el

(setq AnimationEditor
      (ede-cpp-root-project "AnimationEditor" 
							:file "~/Workspace/AnimationEditor/CMakeLists.txt"
                            :include-path'(
										   "SRC/ExternalLib"
										   "SRC/ExternalLib/volumetricMesh-v2.0"
										   "SRC/BaseSim"
										   "SRC/BaseSim/ConAndDrag"
										   "SRC/BaseSim/RedSimulator"
										   "SRC/AnimationEditor"
										   "SRC/AniEditUI"
										   "SRC/IEDS"
										   )

                            :local-variables '( (compile-command . "cd ~/Workspace/AnimationEditor/Build/Release; make -j8 -k") )
                            )
	  )