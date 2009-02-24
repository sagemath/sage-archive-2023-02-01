****************
Calcul distribué
****************

Sage contient une puissante infrastructure de calcul distribué appelé
Distributed Sage (``dsage``).

Vue d'ensemble
==============

Distributed Sage est une infrastructure qui permet de faire du calcul
distribué à partir de l'intérieur de Sage. Il contient un serveur, un client,
des programmes esclaves, et un ensemble de classes que l'on peut dériver
pour écrire des jobs de calcul distribué. Il est destiné en priorité aux
calculs distribués « grossièrement », avec peu de communications entre
les jobs (« grille de calcul »).

Distributed Sage est divisé en trois parties :


#. Le **serveur** est responsable de la distribution des jobs, de leur
   soumission et de leur collecte. Il dispose d'une interface web qui
   permet de surveiller l'état des jobs et d'effectuer diverses autres
   tâches d'aministration.

#. Le **client** a pour rôle de soumettre des jobs au serveur et de
   récupérer les résultats.

#. Les **esclaves** sont les programmes qui font le calcul proprement
   dit.


Prise en main
=============

Voici quelques exemples qui montrent comment démarrer avec ``dsage``.

Exemple 1
---------


#. Lancez ``dsage.setup()``. Cela va configurer la base de données
   SQLite et créer une paire de clés (i.e. une clé publique et une clé
   privée) qui seront utilisées pour la communication SSL. Un
   utilisateur par défaut sera créé, avec comme nom (par défaut) votre
   nom d'utilisateur courant.

#. Lancez ensuite ``d = dsage.start_all()``. Cette commande démarre le
   serveur, le serveur web et :math:`2` esclaves. Elle renvoie un objet
   (``d``) qui représente une connexion vers le serveur. Désormais,
   l'essentiel de vos interactions avec ``dsage`` passeront par cet
   objet ``d``.

#. Ouvrez votre navigateur web à l'adresse http://localhost:8082 pour
   accéder à l'interface web de ``dsage``. À partir de celle-ci, vous
   pourrez voir le statut de vos jobs, les esclaves connectés et
   diverses autres informations importantes sur votre serveur ``dsage``.

#. Commençons par un exemple simple. Tapez ``job = d('2+2')``. En
   consultant l'interface web, vous devriez maintenant voir un nouveau
   job dans le tableau. Un de vos processus esclaves va maintenant
   récupérer le job, l'exécuter et vous renvoyer le résultat. Il est
   possible qu'il ne soit pas encore visible car pour un calcul simple
   comme celui-là, c'est le surcoût correspondant à la communication par
   le réseau qui domine le temps de calcul. Si vous souhaitez attendre
   que votre job se termine, appelez ``job.wait()``, qui bloquera
   jusqu'à ce que le job soit terminé. Vous pourrez ensuite inspecter
   ``job.result`` pour lire le résultat. N'importe quel calcul peut être
   fait de cette façon en appelant ``d``.


Exemple 2
---------

Cet exemple montre comment utiliser la classe ``DistributedFactor``
intégrée à ``dsage``. DistributedFactor tente de factoriser un nombre
par une combinaison de l'algorithme ECM (factorisation par les courbes
elliptiques), du crible quadratique et de l'élimination des petits facteurs
par divisions successives.


#. Lancez ``d = dsage.start_all()`` si vous n'avez pas encore démarré
   votre session ``dsage``. Sinon, vous pouvez continuer à utiliser
   l'instance précédente de ``d``.

#. Démarrez le job de factorisation distribuée par la commande
   ``job_factorisation = DistributedFactor(d, nombre)``. Vous pouvez
   prendre des nombres assez grands, essayez par exemple
   :math:`2^{360}-1`. Pour voir si la factorisation est terminée,
   consultez l'attribut ``job_factorisation.done``. Une fois le job
   terminé, les facteurs premiers trouvés sont disponibles dans
   ``job_factorisation.prime_factors``.


Fichiers
========

``dsage`` garde quelques fichiers dans ``$SAGE_ROOT/.sage/dsage`` :


#. ``pubcert.pem`` et ``cacert.pem`` : la clé publique et la clé privée
   utilisées par le serveur pour les communications SSL.

#. ``dsage_key.pub`` et ``dsage_key`` : les clés utilisées pour
   authentifier l'utilisateur.

#. ``db/`` : sous-répertoire qui contient la base de données ``dsage``.

#. ``*.log`` : journaux du serveur et des esclaves.

#. ``tmp_worker_files/`` : sous-répertoire où les esclaves sauvegardent
   les jobs qu'ils ont traités.

